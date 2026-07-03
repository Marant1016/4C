# Tree-Based Reduced-Lung Solver Implementation Plan

Date: 2026-07-03

## Goal

Implement a custom solver workflow for the reduced-lung model that can eventually use a tree-based linear solver for Newton correction systems.

The intended long-term split is:

| Level | Responsibility |
| --- | --- |
| Nonlinear solver level | Run Newton iterations, evaluate residuals, assemble/update Jacobians, check convergence, update the nonlinear solution. |
| Linear solver level | Solve one linearized Newton system per nonlinear iteration, ideally exploiting the airway-tree structure with Schur complements and bottom-up/top-down elimination. |

The first implementation should avoid changing the lung physics. The safest path is to reuse the existing residual/Jacobian assembly and replace the solver workflow in small steps.

## Current Solver Situation

The current reduced-lung implementation solves a nonlinear algebraic system at each time step:

```text
F(x, time) = 0
```

Here, `x` is the vector of reduced-lung unknowns, mainly element endpoint pressures and flows.

The current workflow is NOX-based:

1. `ReducedLung::reduced_lung_main()` starts the reduced-lung simulation in `src/reduced_lung/src/4C_reduced_lung_main.cpp`.
2. `ReducedLungSimulation::initialize()` builds the topology, model containers, equation ids, maps, vectors, sparse matrix, assembly pipeline, and solver.
3. `ReducedLungSimulation::solve_timestep()` advances time and calls `nox_solver_->solve(current_time_)`.
4. `ReducedLung::NoxSolver` wraps NOX and is defined in `src/reduced_lung/src/4C_reduced_lung_helpers.hpp` and implemented in `src/reduced_lung/src/4C_reduced_lung_helpers.cpp`.
5. `NoxSolver` constructs a `Core::LinAlg::Solver`, binds reduced-lung residual/Jacobian callbacks, and passes everything to `NOX::Nln::Adapter`.
6. NOX manages the nonlinear Newton workflow.
7. Inside each Newton step, NOX calls the 4C linear solver through `NOX::Nln::LinearSystem::apply_jacobian_inverse()`.

NOX currently manages automatically:

| Task | Current owner |
| --- | --- |
| Nonlinear iteration loop | NOX |
| Residual/Jacobian valid-state management | NOX |
| Residual norm status test | NOX |
| Increment norm status test | NOX |
| Maximum nonlinear iteration check | NOX |
| Newton direction computation | NOX |
| Call into `Core::LinAlg::Solver` | NOX linear-system layer |
| Sign convention for Newton direction | NOX solves `J * y = F`, then scales by `-1` |
| Nonlinear solution update | NOX |
| Final solution copy-back | `NOX::Nln::Adapter::solve()` |

Residual assembly currently enters through `ReducedLung::NoxSolver::residual()` and uses the callback pipeline from `create_default_nox_assembly_pipeline()`.

Jacobian assembly currently enters through `ReducedLung::NoxSolver::jacobian()` and uses the same callback pipeline.

The main residual/Jacobian component functions are:

| Component | Residual | Jacobian |
| --- | --- | --- |
| Airways | `Airways::update_residual_vector()` | `Airways::update_jacobian()` |
| Terminal units | `TerminalUnits::update_residual_vector()` | `TerminalUnits::update_jacobian()` |
| Junctions | `Junctions::update_residual_vector()` | `Junctions::update_jacobian()` |
| Boundary conditions | `BoundaryConditions::update_residual_vector()` | `BoundaryConditions::update_jacobian()` |

Important formulation detail: the current reduced-lung model does not use shared nodal unknowns. Each element owns endpoint dofs, and junction equations enforce pressure continuity and flow conservation.

## Proposed Implementation Steps

### Phase 1: Preserve Baseline Behavior

Keep the current NOX path untouched while preparing an alternative path.

1. Use the current NOX workflow as the numerical reference.
2. Reuse existing reduced-lung tests and representative input files as regression targets.
3. Do not change residual equations, Jacobian equations, dof ordering, or row ordering.

### Phase 2: Generalize The Assembly Pipeline

The existing `NoxAssemblyPipeline` is useful outside NOX.

1. Rename or conceptually generalize it to something like `ReducedLungAssemblyPipeline`.
2. Keep the same callback groups: residual assemblers, Jacobian assemblers, state updaters.
3. Keep the current callback order.
4. Keep a compatibility path for `NoxSolver` while the new solver is being developed.

### Phase 3: Add A Custom Newton Solver Beside NOX

Create a reduced-lung Newton solver class parallel to `NoxSolver`.

The custom Newton loop should manually:

1. Set `current_time_`.
2. Synchronize trial `x` into `dofs_` and `locally_relevant_dofs_`.
3. Run state updater callbacks.
4. Assemble residual `F(x)`.
5. Compute residual norm and check `nonlinear_residual_tolerance`.
6. Assemble Jacobian `J(x)`.
7. Solve `J * delta = -F`.
8. Update `x += delta`.
9. Compute increment norm and check `nonlinear_increment_tolerance`.
10. Stop on convergence or throw after `max_nonlinear_iterations`.
11. Synchronize the converged `x` back to reduced-lung dof/state vectors.

First implementation detail: the new Newton solver should call the existing `Core::LinAlg::Solver`. This isolates removal of NOX from development of the custom tree-based linear solver.

### Phase 4: Compare Custom Newton Against NOX

Before writing the tree solver, verify that the custom Newton loop with the existing sparse linear solver reproduces NOX results.

Suggested tests:

1. Single terminal-unit test currently covered by `4C_reduced_lung_nox_solver_test.np2.cpp`.
2. Serial airway with pressure/flow boundary conditions.
3. Bifurcation case with known pressure and flow balances.
4. Terminal-unit cases with both elasticity models if feasible.

Expected agreement should include final solution values, number of time steps completed, and residual/increment convergence behavior. Exact nonlinear iteration counts may differ if the convergence check order differs from NOX.

### Phase 5: Add A Linear Solver Interface For Newton Systems

Introduce a reduced-lung-specific linear solver abstraction used by the custom Newton loop.

The interface should expose a simple operation:

```text
solve(J, F, x, metadata) -> delta
```

Recommended sign convention:

```text
J * delta = -F
x_new = x_old + delta
```

Initial implementations:

1. A wrapper around `Core::LinAlg::Solver` for baseline sparse solves.
2. A custom tree-based solver implementation once metadata and validation are ready.

### Phase 6: Build Tree Metadata

The tree solver needs explicit topology and block layout information, not just a generic sparse matrix.

Build or collect metadata from:

1. Input topology in `ReducedLungParameters::LungTree::Topology`.
2. `first_global_dof_of_ele` and `global_dof_per_ele`.
3. Airway and terminal-unit model data containers.
4. `Junctions::ConnectionData` and `Junctions::BifurcationData`.
5. `BoundaryConditions::BoundaryConditionContainer`.
6. Existing row and dof maps.

Validate at least:

1. Single connected tree if required by the algorithm.
2. One root or clearly defined inlet boundary.
3. Acyclic parent-child topology.
4. Supported branch degree, currently connection or bifurcation.
5. Correct element direction, where `element_nodes[0]` is parent/inlet and `element_nodes[1]` is child/outlet.
6. Boundary conditions sufficient to close the system.

### Phase 7: Implement Tree-Based Linear Solver

Start with a serial implementation and compare against the sparse baseline.

The tree-based linear solver needs input:

| Input | Purpose |
| --- | --- |
| Residual vector `F` or RHS `-F` | Right-hand side of Newton correction system. |
| Jacobian coefficients or structured local derivatives | Linearized element, junction, and boundary equations. |
| Current dof vector `x` | Useful for diagnostics or coefficient recomputation. |
| `dt` and current time | Needed for time-dependent terms and boundary conditions. |
| Dof layout | Map element unknowns to vector entries. |
| Row layout | Map equations to residual/Jacobian rows. |
| Tree topology | Parent-child traversal order. |
| Junction data | Coupling equations between parent and children. |
| Boundary-condition data | Root/leaf constraints. |
| Element model metadata | Block sizes and model types. |

The tree-based linear solver must output:

| Output | Meaning |
| --- | --- |
| Correction vector `delta` | Same global dof numbering as the nonlinear solution vector. |
| Success/failure status | Newton loop must know whether the linear solve succeeded. |
| Optional diagnostics | Linear residual, subtree singularities, pivot information, traversal statistics. |

Implementation strategy:

1. First consume existing assembled sparse `J` plus tree metadata if that is simpler.
2. Compare `delta` against the existing sparse solver on small trees.
3. Then move toward structured block assembly for the Schur complement algorithm.
4. Implement bottom-up elimination from leaves to root.
5. Implement top-down recovery from root to leaves.
6. Add checks for singular or underconstrained subtrees.

### Phase 8: Switch The Reduced-Lung Workflow

After the custom Newton plus tree linear solver path is verified:

1. Add an input option to select NOX, custom Newton with sparse solver, or custom Newton with tree solver.
2. Make the custom path the default only after tests are stable.
3. Remove the NOX dependency from the reduced-lung module only after the NOX path is no longer needed.

## Files, Classes, And Functions Likely To Change

| Area | Likely locations |
| --- | --- |
| Solver construction and per-time-step call | `src/reduced_lung/src/4C_reduced_lung_main.cpp`, especially `ReducedLungSimulation::build_linear_system_and_solver()` and `ReducedLungSimulation::solve_timestep()` |
| Solver wrapper layer | `src/reduced_lung/src/4C_reduced_lung_helpers.hpp` and `src/reduced_lung/src/4C_reduced_lung_helpers.cpp` |
| Assembly pipeline naming/generalization | `NoxAssemblyPipeline`, `NoxSolverContext`, `create_default_nox_assembly_pipeline()` |
| Current NOX wrapper | `ReducedLung::NoxSolver` if kept as one selectable backend |
| New Newton solver | likely new files such as `4C_reduced_lung_newton_solver.hpp/cpp` |
| New reduced-lung linear solver interface | likely new files such as `4C_reduced_lung_linear_solver.hpp/cpp` |
| New tree metadata builder | likely new files such as `4C_reduced_lung_tree_structure.hpp/cpp` |
| New tree linear solver | likely new files such as `4C_reduced_lung_tree_linear_solver.hpp/cpp` |
| Input options | `src/reduced_lung/src/4C_reduced_lung_input.hpp` and `src/reduced_lung/src/4C_reduced_lung_input.cpp` |
| Build dependencies | `src/reduced_lung/CMakeLists.txt`, especially removal of `solver_nonlin_nox` after the NOX path is removed |
| Tests | `src/reduced_lung/tests/`, especially adding tests beside `4C_reduced_lung_nox_solver_test.np2.cpp` |

Important current source locations:

| Topic | Location |
| --- | --- |
| Main reduced-lung entry point | `ReducedLung::reduced_lung_main()` in `4C_reduced_lung_main.cpp` |
| Time-step solve call | `ReducedLungSimulation::solve_timestep()` in `4C_reduced_lung_main.cpp` |
| Vector/matrix/solver setup | `ReducedLungSimulation::build_linear_system_and_solver()` in `4C_reduced_lung_main.cpp` |
| NOX wrapper | `ReducedLung::NoxSolver` in `4C_reduced_lung_helpers.hpp/cpp` |
| NOX parameter list | `create_nox_parameter_list()` in `4C_reduced_lung_helpers.cpp` |
| Residual callback entry | `NoxSolver::residual()` in `4C_reduced_lung_helpers.cpp` |
| Jacobian callback entry | `NoxSolver::jacobian()` in `4C_reduced_lung_helpers.cpp` |
| State synchronization | `NoxSolver::sync_state_from_x()` in `4C_reduced_lung_helpers.cpp` |
| Current linear solver call | `NOX::Nln::LinearSystem::apply_jacobian_inverse()` in `src/solver_nonlin_nox/4C_solver_nonlin_nox_linearsystem.cpp` |

## Files, Classes, And Functions To Avoid Changing At First

Avoid changing these during the first solver refactor:

| Area | Reason |
| --- | --- |
| Airway residual/Jacobian physics | Already separated from NOX and tested. |
| Terminal-unit residual/Jacobian physics | Already separated from NOX and tested. |
| Junction residual/Jacobian equations | Needed as the baseline tree coupling equations. |
| Boundary-condition residual/Jacobian equations | Needed to close the system and compare against NOX. |
| Dof ordering | Many helper functions assume fixed element-local offsets. |
| Row ordering | Residual/Jacobian assembly assumes current equation ordering. |
| Global problem-type switch | Not needed for replacing the reduced-lung solver internals. |
| `solver_nonlin_nox` internals | Avoid changing global NOX infrastructure for a reduced-lung-specific solver. |
| `src/reduced_lung/src/1d_pipe_flow/*` | Separate reduced-lung 1D pipe-flow implementation. |
| Model registries | Not needed unless adding new physics models. |

Current model physics files to leave alone initially:

| Component | Files |
| --- | --- |
| Airways | `airways/4C_reduced_lung_airways.cpp`, `airways/4C_reduced_lung_airways_wall_mechanics.cpp`, `airways/4C_reduced_lung_airways_flow_resistance.cpp` |
| Terminal units | `terminal_units/4C_reduced_lung_terminal_unit.cpp`, `terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp`, `terminal_units/4C_reduced_lung_terminal_unit_elasticity.cpp` |
| Junctions | `4C_reduced_lung_junctions.cpp` |
| Boundary conditions | `4C_reduced_lung_boundary_conditions.cpp` |

## Risks And Open Questions

Risks:

1. The current `x_` vector is built on `row_map_`, while `dofs_` is built on `locally_owned_dof_map_`. The current implementation assumes compatible ordering for `Core::LinAlg::export_to(x, dofs_)`.
2. Several Jacobian assemblers rely on sparse-matrix fill state. First assembly inserts entries; later assemblies replace nonlinear entries or skip constant blocks.
3. Junction creation depends on element direction. Wrong `element_nodes` ordering can produce wrong parent/child interpretation.
4. Boundary-condition ownership is parallel-distribution dependent.
5. A tree solver may need a stricter topology definition than the current sparse solver.
6. The current duplicated endpoint dof formulation may not be the cleanest formulation for a Schur complement tree solver.
7. Removing NOX also removes NOX convergence behavior, failure handling, and line-search infrastructure.
8. Full-step Newton may be less robust than NOX if damping or line search becomes necessary.
9. Parallel tree elimination may require nontrivial communication if subtrees are distributed over MPI ranks.

Open questions to discuss before implementation:

1. Should the thesis solver preserve the current duplicated endpoint dof formulation or move to shared node variables?
2. Is serial-only tree solving acceptable for the first implementation?
3. Should unsupported topologies fall back to the sparse solver or throw a clear error?
4. Is full-step Newton acceptable, or is damping/line search required?
5. Which model variants must the tree solver support first?
6. Should the first tree solver consume assembled sparse `J`, or should it assemble structured element/junction blocks directly?
7. What is the primary validation target: agreement with NOX/UMFPACK, performance, or physiological benchmark behavior?
8. Should solver selection be an input-file option from the beginning?
9. How strict should the tree topology validation be for real lung trees?

## Checklist For Future Sessions

Before making code changes:

1. Read `src/reduced_lung/AGENTS.md`.
2. Read `src/reduced_lung/info_oc/reduced_lung_baseline.md`.
3. Read this file.
4. Confirm whether the task is documentation-only or implementation.
5. Avoid touching `src/reduced_lung/src/1d_pipe_flow/*` unless explicitly requested.
6. Check current worktree status and avoid reverting user changes.

When implementing:

1. Keep the NOX path working until the custom path is verified.
2. Reuse the existing residual/Jacobian assembly first.
3. Add the custom Newton loop before adding the tree linear solver.
4. Use the existing sparse linear solver as the first backend for the custom Newton loop.
5. Compare custom Newton results against NOX on small tests.
6. Introduce tree metadata and validation before implementing Schur complement elimination.
7. Compare tree solver corrections against the sparse solver on small systems.
8. Only remove the NOX dependency after the replacement path is stable.

Definition of a safe first milestone:

1. A selectable custom Newton solver exists beside `NoxSolver`.
2. It reuses the existing assembly callbacks.
3. It calls `Core::LinAlg::Solver` for `J * delta = -F`.
4. It reproduces at least the single-terminal-unit NOX test.
5. No airway, terminal-unit, junction, or boundary physics code has been changed.
