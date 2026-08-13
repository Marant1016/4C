# reduced_lung OpenCode Instructions

## Scope

- These instructions apply to `src/reduced_lung/` and supplement the repository-level `AGENTS.md` in the 4C source root.
- Treat the current source as authoritative. Historical implementation notes may exist locally, but they can be stale after later solver and cleanup changes.
- Keep `src/1d_pipe_flow/` separate. It is a different reduced-lung 1D pipe-flow path; do not edit it unless the task explicitly names it.

## Module Context

- `reduced_lung` implements a reduced-dimensional lung-tree solver for `PROBLEMTYPE: "Reduced_Lung"` inputs.
- The model uses a duplicated-endpoint formulation: each element owns its endpoint pressure/flow dofs, and junction equations enforce pressure continuity and flow conservation between elements.
- Elements are either `Airway` or `TerminalUnit`. Airways support rigid and Kelvin-Voigt wall mechanics plus linear and nonlinear flow resistance. Terminal units support Kelvin-Voigt and four-element Maxwell rheology plus linear and Ogden elasticity.
- Input topology is directional. In `element_nodes`, the first node is the inlet/parent side and the second node is the outlet/child side. This ordering is critical for junction creation and the tree solver.

## Current Solver Policy

- Runtime solver selection is controlled by `reduced_dimensional_lung.dynamics.nonlinear_solver` in the input file.
- Supported values are `Nox`, `NewtonSparse`, and `NewtonTree`.
- `Nox` is the default/reference workflow and remains MPI-capable through the existing NOX infrastructure.
- `NewtonSparse` is the custom full-step Newton workflow using `SparseNewtonLinearSolver` and the regular 4C sparse linear solver backend.
- `NewtonTree` is the optimized custom full-step Newton workflow using the serial `TreeNewtonLinearSolver` with structured tree-block linearization.
- `NewtonTree` is serial-only. Running it with more than one MPI rank should fail fast; do not add silent fallback to `Nox` or `NewtonSparse`, because that would make performance runs misleading.

## Important Files

- `src/4C_reduced_lung_main.cpp`: top-level simulation setup, time loop, solver selection, serial-only `NewtonTree` check, and profile summary printing.
- `src/4C_reduced_lung_input.hpp` and `src/4C_reduced_lung_input.cpp`: `ReducedLungParameters`, valid input parameters, and the `NonlinearSolverType` enum.
- `src/4C_reduced_lung_helpers.hpp` and `src/4C_reduced_lung_helpers.cpp`: discretization construction, dof/row/column maps, model creation helpers, NOX wrapper, output collection, and `ReducedLungAssemblyPipeline` registration.
- `src/solver/4C_reduced_lung_newton_solver.hpp` and `src/solver/4C_reduced_lung_newton_solver.cpp`: custom full-step Newton loop used by `NewtonSparse` and `NewtonTree`.
- `src/solver/4C_reduced_lung_newton_linear_solver.hpp` and `src/solver/4C_reduced_lung_newton_linear_solver.cpp`: Newton linear-solver interface and sparse backend.
- `src/solver/tree/4C_reduced_lung_tree_metadata.hpp` and `src/solver/tree/4C_reduced_lung_tree_metadata.cpp`: tree metadata construction and validation.
- `src/solver/tree/4C_reduced_lung_tree_linearization.hpp` and `src/solver/tree/4C_reduced_lung_tree_linearization.cpp`: structured coefficient target and row-oriented tree linearization storage.
- `src/solver/tree/4C_reduced_lung_tree_linear_solver.hpp` and `src/solver/tree/4C_reduced_lung_tree_linear_solver.cpp`: optimized serial bottom-up/top-down tree Newton linear solver.
- `src/solver/4C_reduced_lung_solver_profiles.hpp`: profiling counters for custom Newton, sparse Newton, and tree linear solve phases.
- `src/airways/`, `src/terminal_units/`, `src/4C_reduced_lung_junctions.*`, and `src/4C_reduced_lung_boundary_conditions.*`: physics residuals, sparse Jacobians, structured tree coefficients, state updates, and output callbacks.

## Setup Workflow

- `ReducedLung::reduced_lung_main(Global::Problem&)` creates a `ReducedLungSimulation`, initializes it, and runs the time loop.
- Initialization validates dynamics parameters, builds the lightweight line-element discretization from input topology, creates element model containers, builds node-level boundary conditions and junctions, assigns equation ids, creates maps/local ids, and then builds the selected nonlinear solver.
- Dof ordering is fixed by element type: rigid airway `p1, p2, q1`; Kelvin-Voigt airway `p1, p2, q1, q2`; terminal unit `p1, p2, q`.
- Row ordering is state equations first, then junction/coupling equations, then boundary-condition equations.
- The solution vector `x_` is created on `row_map_`, while owned dofs use `locally_owned_dof_map_`. Existing code assumes the equation and unknown layouts are compatible for synchronization. Be careful before changing row or dof ordering.

## Assembly Pipeline

- `create_default_reduced_lung_assembly_pipeline()` registers residual, sparse Jacobian, structured tree-linearization, static tree-linearization, tree-capacity, and state-update callbacks.
- Residual callback order is airways, terminal units, junctions, then boundary conditions.
- Sparse Jacobian callback order is airways, terminal units, junctions, then boundary conditions.
- Structured tree coefficients keep physics derivatives in the physics modules. Do not duplicate airway, terminal-unit, junction, or boundary derivative formulas inside the tree solver.
- Static tree coefficients are appended once when possible; dynamic coefficients are replaced through the `TreeCoefficientAssemblyTarget` interface.
- `TreeNewtonLinearSolver::direct_tree_coefficient_target()` lets the runtime `NewtonTree` path write directly into solver coefficient arrays instead of rebuilding a generic `TreeLinearization`. Keep this fast path intact unless a task is explicitly about changing structured assembly.

## Newton Workflows

- `NoxSolver` wraps the existing 4C NOX adapter and uses `Core::LinAlg::Solver` through NOX.
- `NewtonSolver` owns the custom full-step Newton loop. It synchronizes `x` to owned and locally relevant dofs, updates model internal state, assembles residuals, checks residual convergence, assembles either sparse Jacobian or structured tree coefficients, solves `J * delta = -F`, and applies `x += delta`.
- The residual norm is the authoritative convergence check. If a nonzero iteration has a small increment but the residual is not converged, `NewtonSolver` treats that as stagnation.
- For `NewtonSparse`, `NewtonSolver` assembles the sparse Jacobian and calls `SparseNewtonLinearSolver`.
- For `NewtonTree`, `NewtonSolver` assembles structured tree coefficients and calls `TreeNewtonLinearSolver`; the sparse Jacobian should not be assembled in the runtime tree path.

## Tree Solver Context

- The tree solver preserves the same reduced-lung equations and solves the same Newton correction system as the sparse solver: `J * delta = -F`, followed by `x += delta`.
- The serial tree algorithm condenses subtrees bottom-up into affine inlet relations `delta_q_in = G * delta_p_in + h`, closes the root with the root inlet pressure boundary row, then recovers all dof corrections top-down.
- Tree metadata validates assumptions required by this algorithm: single rooted directed tree, no cycles, all elements reachable from the root, branch degree at most two, matching connection/bifurcation metadata, terminal units as leaves, square system, root inlet boundary, and outlet boundary for each leaf.
- The current root closure effectively expects a root inlet pressure boundary. A root flow closure would require additional solver logic.
- Current production element block sizes are `2x2` for rigid airways and terminal units, and `3x3` for Kelvin-Voigt airways. Optimized batch/SIMD-style paths exist for these common shapes, with generic dense fallback for unsupported or numerically difficult cases.
- `TreeNewtonLinearSolverCoefficientSource::SparseJacobian` remains useful for validation, but runtime `NewtonTree` should use `StructuredTreeBlocks`.

## Performance And Profiling

- Enable reduced-lung tree profiling with `FOUR_C_REDUCED_LUNG_TREE_PROFILE=1`.
- Profile output is printed from `ReducedLungSimulation::print_tree_profile_summary()` and includes residual phase timings, structured tree assembly timings, tree solve timings, dense solve counts, SIMD/group counters, fallback counts, element counts, and max block size.
- For large serial `NewtonTree` cases, pay attention to `residual_s`, `tree_assembly_s`, `tree_solve_s`, `tree_dense_fallbacks`, `tree_unsupported_fallbacks`, `tree_simd_lanes`, and `tree_max_block`.
- Normal large-tree production runs should avoid sparse Jacobian assembly, unsupported block fallbacks, and dense fallbacks unless the case genuinely needs a numerical fallback.

## Inputs And Archived Fixtures

- Reduced-lung runtime and benchmark inputs created during the custom Newton work are archived outside the source tree under `/scratch/Rodriguez/workspace/4C/files/`.
- Many archived YAMLs use relative `from_file` JSON references. Run them from the directory containing the YAML file.
- Representative archived unit inputs are in `/scratch/Rodriguez/workspace/4C/files/unit_tests/`, including small `NewtonTree`, `NewtonSparse`, and `Nox` variants plus a gen10 tree case.
- Generated gen16 benchmark inputs are in `/scratch/Rodriguez/workspace/4C/files/gen16_inputs/`.
- Benchmark scripts are in `/scratch/Rodriguez/workspace/4C/files/benchmarks/scripts/`.
- Archived GoogleTest sources for custom Newton/tree work are under `/scratch/Rodriguez/workspace/4C/files/google_tests/reduced_lung/`. They are not auto-discovered from that location; restore them into the test tree only if a task explicitly asks for that.

## Verification

- From the 4C source root, prefer targeted reduced-lung builds and tests when a local build exists.
- Useful build commands are `cmake --build build/debug --target reduced_lung_objs --parallel 4` and `cmake --build build/debug --target unittests_reduced_lung unittests_reduced_lung.np2 --parallel 4`.
- Useful CTest commands are `ctest -R "^unittests_reduced_lung(\.np2)?$" --output-on-failure` and targeted reduced-lung input regexes such as `ctest -R "reduced_lung_.*newton_.*\.4C\.yaml-p1$" --output-on-failure` when those inputs are registered in the local checkout.
- To run archived gen16 inputs directly, use a release 4C executable when performance matters, run from `/scratch/Rodriguez/workspace/4C/files/gen16_inputs/`, and keep `NewtonTree` serial.
- Example direct run: `FOUR_C_REDUCED_LUNG_TREE_PROFILE=1 /scratch/Rodriguez/workspace/4C/4C/build/release/4C reduced_lung_lung_tree_gen16_500steps_newton_tree.4C.yaml /scratch/Rodriguez/workspace/4C/output/gen16_newton_tree_profile`.

## Working Rules For Future Changes

- Keep `Nox` working as the safe default/reference path unless the task explicitly removes it.
- Keep `NewtonSparse` as the custom Newton sparse reference path for comparing against `NewtonTree`.
- Keep `NewtonTree` explicit opt-in and serial-only.
- Do not change residual equations, Jacobian formulas, dof ordering, row ordering, or topology directionality as part of performance-only work.
- When optimizing, preserve a validation path that can compare tree corrections against sparse corrections on small cases.
- Be careful with sparse matrix fill state: some Jacobian assemblers insert on the first unfilled pass and replace or skip entries after `SparseMatrix::complete()`.
- Be careful with model containers: airways and terminal units are grouped by model variants, and per-element arrays inside each group must remain aligned.
- Add or update tests when behavior changes. For performance-only refactors, verify numerical equivalence against the existing sparse or NOX reference path whenever practical.
