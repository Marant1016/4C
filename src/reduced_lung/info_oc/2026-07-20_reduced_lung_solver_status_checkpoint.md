# Reduced-Lung Solver Status Checkpoint - 2026-07-20

This checkpoint summarizes the current reduced-lung solver implementation after completing the phased tree-solver workflow plan. It is written for readers who know the original NOX-based reduced-lung implementation and need to understand what changed, what now runs, and what remains future work.

The important high-level point is:

```text
The reduced-lung module now has selectable nonlinear solver workflows:
Nox, NewtonSparse, and NewtonTree.
```

NOX remains available and remains the default. The new custom Newton and tree linear solver paths are implemented beside it.

## Original NOX Baseline

Originally, the reduced-lung time-step workflow was entirely NOX-based.

The runtime path was:

```text
ReducedLung::reduced_lung_main()
  -> ReducedLungSimulation::initialize()
  -> ReducedLungSimulation::run()
  -> ReducedLungSimulation::solve_timestep()
  -> NoxSolver::solve(current_time)
```

`NoxSolver` wrapped the global NOX nonlinear solver infrastructure. It owned or referenced the solution vector, reduced-lung dof vectors, assembly callbacks, sparse Jacobian operator, and a `Core::LinAlg::Solver` instance.

NOX handled:

- nonlinear iteration control,
- residual and Jacobian callback invocation,
- convergence checks,
- Newton direction computation,
- calls into the 4C sparse linear solver,
- final solution copy-back.

The reduced-lung physics assembly was already separated into component callbacks:

- airways residual/Jacobian assembly,
- terminal-unit residual/Jacobian assembly,
- junction residual/Jacobian assembly,
- boundary-condition residual/Jacobian assembly.

That original assembly behavior is still preserved.

## Runtime Solver Options

The reduced-lung input now has a nonlinear solver workflow selector:

```yaml
reduced_dimensional_lung:
  dynamics:
    nonlinear_solver: Nox
```

The corresponding C++ enum is:

```text
ReducedLungParameters::NonlinearSolverType
```

Supported values are:

| Option | Nonlinear Solver | Linear Correction Solver | Current Role |
| --- | --- | --- | --- |
| `Nox` | NOX | `Core::LinAlg::Solver` through NOX | Default production/reference path |
| `NewtonSparse` | `NewtonSolver` | `SparseNewtonLinearSolver` | Custom Newton validation path using sparse backend |
| `NewtonTree` | `NewtonSolver` | `TreeNewtonLinearSolver` | New tree-based linear solver path |

The default is:

```text
Nox
```

Existing input files that omit `nonlinear_solver` therefore keep the original NOX behavior.

## Runtime Construction

Solver selection happens in `ReducedLungSimulation::build_linear_system_and_solver()` in:

```text
src/reduced_lung/src/4C_reduced_lung_main.cpp
```

The setup still builds the common reduced-lung objects first:

- discretization,
- airway and terminal-unit model containers,
- junction and boundary-condition containers,
- equation ids,
- dof maps,
- row and column maps,
- solution/residual vectors,
- sparse Jacobian matrix,
- assembly pipeline.

After that, it selects the nonlinear workflow:

```text
Nox
  -> build_nox_solver()

NewtonSparse
  -> build_newton_solver_with_sparse_linear_solver()
  -> build_newton_solver()

NewtonTree
  -> build_newton_solver_with_tree_linear_solver()
  -> build_reduced_lung_tree_metadata(...)
  -> build_newton_solver()
```

`ReducedLungSimulation::solve_timestep()` then calls either:

```text
NoxSolver::solve(current_time)
```

or:

```text
NewtonSolver::solve(current_time)
```

The end-of-step model updates are unchanged:

```text
TerminalUnits::end_of_timestep_routine(...)
Airways::end_of_timestep_routine(...)
```

## Main Files And Classes Introduced

### Assembly Pipeline

Files:

```text
src/reduced_lung/src/4C_reduced_lung_helpers.hpp
src/reduced_lung/src/4C_reduced_lung_helpers.cpp
```

Main type:

```text
ReducedLungAssemblyPipeline
```

This generalized the previous NOX-specific assembly callback grouping. It stores:

- residual assemblers,
- Jacobian assemblers,
- state updaters.

The compatibility alias remains:

```text
NoxAssemblyPipeline = ReducedLungAssemblyPipeline
```

### Custom Newton Solver

Files:

```text
src/reduced_lung/src/4C_reduced_lung_newton_solver.hpp
src/reduced_lung/src/4C_reduced_lung_newton_solver.cpp
```

Main types:

```text
NewtonSolverContext
NewtonSolver
```

`NewtonSolver` owns the custom full-step Newton loop and uses a `NewtonLinearSolver` for each correction system.

### Newton Linear Solver Interface

Files:

```text
src/reduced_lung/src/4C_reduced_lung_linear_solver.hpp
src/reduced_lung/src/4C_reduced_lung_linear_solver.cpp
```

Main types:

```text
NewtonLinearSystemMetadata
NewtonLinearSolver
SparseNewtonLinearSolverContext
SparseNewtonLinearSolver
```

`NewtonLinearSolver` is the abstraction used by `NewtonSolver`.

`SparseNewtonLinearSolver` wraps `Core::LinAlg::Solver` and preserves the sparse backend behavior.

### Tree Metadata

Files:

```text
src/reduced_lung/src/4C_reduced_lung_tree_metadata.hpp
src/reduced_lung/src/4C_reduced_lung_tree_metadata.cpp
```

Main types:

```text
TreeElementMetadata
TreeJunctionMetadata
TreeBoundaryConditionMetadata
ReducedLungTreeMetadata
ReducedLungTreeMetadataContext
```

Main builder:

```text
build_reduced_lung_tree_metadata(context)
```

The metadata builder creates the explicit tree and layout data needed by the tree linear solver.

### Tree Linear Solver

Files:

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
```

Main types:

```text
TreeNewtonLinearSolverContext
TreeNewtonLinearSolver
```

`TreeNewtonLinearSolver` implements `NewtonLinearSolver` and performs the serial bottom-up/top-down tree correction solve.

### Input And Runtime Switch

Files:

```text
src/reduced_lung/src/4C_reduced_lung_input.hpp
src/reduced_lung/src/4C_reduced_lung_input.cpp
src/reduced_lung/src/4C_reduced_lung_main.cpp
```

Main addition:

```text
ReducedLungParameters::NonlinearSolverType nonlinear_solver
```

This selects `Nox`, `NewtonSparse`, or `NewtonTree` at runtime.

## Custom Newton Algorithm

The custom Newton algorithm is implemented in `NewtonSolver::solve()`.

The mathematical convention is:

```text
J * delta = -F
x_new = x_old + delta
```

For each physical time step, `ReducedLungSimulation::solve_timestep()` advances time and calls `NewtonSolver::solve(current_time)` when a custom Newton workflow is selected.

Inside the nonlinear iteration loop, `NewtonSolver` does the following:

1. Synchronize the current nonlinear solution vector `x_solution_` into the owned dof vector `dofs_`.
2. Export `dofs_` into `locally_relevant_dofs_`.
3. Run all registered state updater callbacks in the assembly pipeline.
4. Zero and assemble the residual vector using all residual assembler callbacks.
5. Compute the residual 2-norm.
6. Check residual convergence against `nonlinear_residual_tolerance`.
7. Check increment convergence against `nonlinear_increment_tolerance`; the first iteration skips the increment check by treating it as converged.
8. If both residual and increment tests pass, return the number of Newton corrections applied.
9. If the maximum nonlinear iteration count has been reached, throw a non-convergence error.
10. Assemble the Jacobian using all Jacobian assembler callbacks.
11. Complete the sparse matrix if it has not been completed yet.
12. Build `NewtonLinearSystemMetadata` with current time, time-step size, and nonlinear iteration index.
13. Call the selected `NewtonLinearSolver` to compute `delta`.
14. Compute the increment 2-norm from `delta`.
15. Update the nonlinear solution with `x_solution_.update(1.0, delta_, 1.0)`.

This is a full-step Newton method. There is currently no damping or line search in `NewtonSolver`.

## Sparse Newton Linear Solver

`SparseNewtonLinearSolver` is the sparse backend for the custom Newton loop.

It receives:

- completed or completable sparse Jacobian,
- residual vector,
- current nonlinear solution vector,
- Newton metadata,
- correction vector `delta`.

It forms:

```text
rhs = -residual
```

Then it calls:

```text
Core::LinAlg::Solver::solve(jacobian, delta, rhs, solver_params)
```

It sets:

```text
refactor = true
reset = metadata.nonlinear_iteration == 0
```

and forwards a configured projector if present in the linear solver parameters.

This path is useful because it isolates custom Newton behavior from tree-linear-solver behavior.

## Tree Metadata

`ReducedLungTreeMetadata` is built from existing reduced-lung setup objects, not from a new input format.

It uses:

- input topology from `ReducedLungParameters::LungTree::Topology`,
- `first_global_dof_of_ele`,
- `global_dof_per_ele`,
- airway model data containers,
- terminal-unit model data containers,
- `Junctions::ConnectionData`,
- `Junctions::BifurcationData`,
- `BoundaryConditions::BoundaryConditionContainer`,
- row map,
- locally relevant dof map.

For each tree element, metadata records:

- global element id,
- element kind, airway or terminal unit,
- inlet and outlet node ids,
- parent element index,
- up to two child element indices,
- first global dof,
- number of dofs,
- global dof ids,
- local dof ids,
- first local/global state-equation ids,
- number of state equations.

For junctions, it records:

- connection or bifurcation kind,
- parent element index,
- child element indices,
- first local/global equation ids,
- number of equations,
- global/local dof ids involved in the junction equations.

For boundary conditions, it records:

- pressure or flow type,
- inlet or outlet side,
- node id,
- element index,
- local/global equation id,
- global/local constrained dof id.

It also records:

- root element index,
- root node id,
- airway element indices,
- terminal-unit element indices,
- bottom-up traversal layers,
- top-down traversal layers,
- global dof count,
- global equation count,
- locally relevant dof count.

The metadata builder validates the directed tree assumptions before the tree solver is constructed. Current validation includes:

- one connected tree,
- exactly one root,
- acyclic directed parent-child topology,
- branch degree at most two,
- terminal units are leaves,
- connection and bifurcation metadata match directed topology,
- every parent with children has matching junction metadata,
- root inlet has a boundary condition,
- every leaf outlet has a boundary condition,
- equation count equals dof count.

## TreeNewtonLinearSolver: Implemented Behavior

`TreeNewtonLinearSolver` currently solves one Newton correction system for the custom Newton loop.

It implements:

```text
NewtonLinearSolver::solve(jacobian, residual, x, metadata, delta)
```

The solver is serial-only and requires:

- a completed sparse Jacobian,
- the residual vector,
- `ReducedLungTreeMetadata`,
- all residual rows locally available,
- all correction dofs locally available,
- a root inlet boundary condition,
- one outlet boundary condition for each leaf element.

The current implementation still consumes the assembled sparse Jacobian and residual. It does not call UMFPACK or `Core::LinAlg::Solver`. Instead, it extracts the local matrix coefficients it needs from the sparse matrix and applies a tree-structured solve.

### Coefficient Source

The current tree solver reads coefficients from the assembled sparse Jacobian with row/column lookups.

This means the path is currently:

```text
existing residual/Jacobian assembly
  -> sparse J and residual F
  -> TreeNewtonLinearSolver extracts coefficients from J
  -> tree bottom-up/top-down solve
  -> delta
```

This is an intermediate, validation-oriented design. It avoids duplicating all airway, terminal-unit, junction, and boundary derivative formulas inside the tree solver.

### Bottom-Up Subtree Condensation

The solver processes elements in `bottom_up_layers`, from leaves toward the root.

For every element, it constructs a local affine relation between inlet pressure correction and inlet flow correction:

```text
delta_q_in = G * delta_p_in + h
```

This relation is stored as a `SubtreeRelation`.

The element-local unknowns are all element dofs except the inlet pressure dof. The inlet pressure correction is treated as the remaining external parameter of the subtree.

For a leaf element, the local system uses:

- the element state-equation rows,
- the leaf outlet boundary-condition row.

For an internal element, the local system uses:

- the element state-equation rows,
- the downstream junction flow-conservation row.

For an internal element, each child subtree has already been condensed into its own relation:

```text
child_delta_q_in = child_G * child_delta_p_in + child_h
```

The pressure-continuity junction rows express child inlet pressure corrections in terms of the parent outlet pressure correction. The child flow relation is then substituted into the parent flow-conservation row.

### Schur-Complement Interpretation

The bottom-up pass is a Schur-complement elimination over the directed tree.

At each subtree:

- child internal dofs are already eliminated,
- child behavior is represented only by a scalar inlet relation,
- the parent eliminates its own local dofs except inlet pressure,
- the result is a new scalar relation for the parent subtree.

This is the C++ implementation of the same core idea used in the notebook prototype:

```text
compress each downstream subtree to dQ = G dP + h
```

The current C++ implementation generalizes that idea by extracting coefficients from the actual assembled 4C reduced-lung Jacobian rather than hard-coding one model's derivative formulas.

### Local Dense Solves

Each element condensation forms a small dense local system.

The solver solves that local dense system twice:
```text
constant/intercept solve -> h contribution
slope solve              -> G contribution
```

The current implementation uses a small internal Gaussian-elimination helper with partial pivoting.

If a local pivot is smaller than the configured tolerance, the solver throws an error for a singular or underconstrained local block.

### Root System Solution

After bottom-up condensation reaches the root, the root inlet boundary condition determines the root inlet pressure correction.

In the current implementation, this is a scalar equation from the root inlet boundary row:

```text
boundary_coefficient * delta_p_root = -F_boundary
```

So:

```text
delta_p_root = -F_boundary / boundary_coefficient
```

This root value starts the top-down recovery.

### Top-Down Correction Recovery

The solver processes elements in `top_down_layers`, from root to leaves.

Each element stores `ElementRecoveryData` from the bottom-up phase. This recovery data expresses all element-local unknown corrections as affine functions of that element's inlet pressure correction.

For each element:

1. The inlet pressure correction is known from either the root boundary or the parent junction recovery.
2. The solver recovers all element-local dof corrections.
3. The recovered values are written into the global correction vector `delta` using global dof ids.
4. If the element has children, pressure-continuity rows are used to compute each child inlet pressure correction.

### Global Correction Vector

The final output is the global Newton correction vector:

```text
delta
```

It uses the same global dof numbering and map convention as the custom Newton sparse path.

After `TreeNewtonLinearSolver::solve()` returns, `NewtonSolver` applies:

```text
x_solution_ = x_solution_ + delta
```

## Verified And Working

The implementation has these verification points documented in the phase notes:

- custom Newton analytical single terminal-unit test,
- custom Newton versus NOX comparisons,
- sparse Newton linear solver wrapper tests through custom Newton workflows,
- tree metadata construction tests,
- tree linear solver correction comparisons against the sparse solver,
- runtime option default test,
- representative reduced-lung input tests using the default NOX path.

Important test files include:

```text
src/reduced_lung/tests/4C_reduced_lung_newton_solver_test.np2.cpp
src/reduced_lung/tests/4C_reduced_lung_newton_vs_nox_test.np2.cpp
src/reduced_lung/tests/4C_reduced_lung_tree_metadata_test.cpp
src/reduced_lung/tests/4C_reduced_lung_tree_linear_solver_test.cpp
src/reduced_lung/tests/4C_reduced_lung_input_pipeline_test.cpp
```

The standard reduced-lung unit-test command covering the new unit tests is:

```text
ctest -R "^unittests_reduced_lung(\.np2)?$" --output-on-failure
```

This runs:

```text
unittests_reduced_lung
unittests_reduced_lung.np2
```

The serial target includes the tree metadata, tree linear solver, and input default tests. The `.np2` target includes the two-process Newton/NOX tests.

The Phase 8 verification also ran representative existing runtime input tests that still use the default NOX path:

```text
reduced_lung_terminal_unit.4C.yaml-p1
reduced_lung_serial_airways_flow.4C.yaml-p1
reduced_lung_aw_bifurcation_flow.4C.yaml-p1
```

## Current Limitations

The new workflow is implemented, but the final high-performance tree solver is not complete yet.

Implemented behavior:

- `Nox`, `NewtonSparse`, and `NewtonTree` are selectable runtime workflows.
- `Nox` remains the default.
- `NewtonSparse` uses custom Newton with the 4C sparse solver.
- `NewtonTree` uses custom Newton with the serial tree linear solver.
- `NewtonTree` performs bottom-up/top-down tree correction solves.

Current limitations:

- `NewtonTree` is serial-only.
- `NewtonTree` still requires the assembled sparse Jacobian as its coefficient source.
- The tree solver does not yet consume structured derivative blocks directly.
- The tree solver is not yet the default.
- NOX remains a required dependency because it is still selectable and still the default.
- There is no damping or line search in the custom Newton solver.
- The tree solver has been directly compared against the sparse solver on selected small cases, not every supported model combination.
- Explicit tree-solver validation is still needed for Kelvin-Voigt airway cases, nonlinear airway resistance, Four-element Maxwell terminal units, Ogden terminal units in more topologies, and larger realistic trees.
- Parallel tree elimination and communication are not implemented.
- The current duplicated endpoint dof formulation is preserved; no shared-node reformulation has been attempted.

## Recommended Next Steps

The next work should separate two goals: robustness validation and performance optimization.

### Validate The Existing Tree Path

Recommended validation work:

- Add tests for `NewtonTree` runtime input files, not only direct linear correction tests.
- Compare complete time-step results for `NewtonSparse` and `NewtonTree`.
- Compare `Nox`, `NewtonSparse`, and `NewtonTree` on all representative reduced-lung YAML inputs where serial execution is acceptable.
- Add explicit tree-solver tests for Kelvin-Voigt airways.
- Add explicit tree-solver tests for nonlinear airway resistance.
- Add explicit tree-solver tests for Four-element Maxwell terminal units.
- Add larger tree tests with mixed airways and terminal units.
- Decide whether unsupported tree layouts should throw or fall back to `NewtonSparse`.

### Move Toward Structured Local Blocks

For performance, the main improvement is to avoid building the generic global sparse Jacobian for the `NewtonTree` path.

The target should become:

```text
current target path:
  residual assembly
  sparse Jacobian assembly
  sparse coefficient extraction
  tree solve

future high-performance path:
  residual assembly
  structured local derivative/block extraction
  tree solve
```

Recommended design direction:

- Keep physics derivatives in airway, terminal-unit, junction, and boundary-condition modules.
- Add structured derivative block callbacks beside the sparse Jacobian callbacks.
- Represent element, junction, and boundary derivatives as small dense blocks tied to tree metadata.
- Let `TreeNewtonLinearSolver` consume these blocks directly.
- Keep the sparse Jacobian path as a reference and fallback.

### Reduce Allocations And Reuse Symbolic Information

The current tree solver creates work vectors and dense local systems inside each solve.

Performance-oriented improvements:

- Preallocate one work object per element.
- Store local row/dof lookup data once in tree metadata or a derived symbolic tree-solver plan.
- Reuse local dense matrix storage across Newton iterations.
- Reuse symbolic child/parent coupling patterns.
- Separate symbolic preprocessing from numeric coefficient updates.
- Avoid repeated sparse row searches once structured blocks exist.

### Benchmark And Profile

Recommended measurement work:

- Benchmark `Nox`, `NewtonSparse`, and `NewtonTree` on the same serial inputs.
- Measure assembly time separately from linear solve time.
- Measure sparse matrix assembly/completion cost separately from coefficient extraction and tree solve cost.
- Profile allocations in the tree solver.
- Profile local dense solve cost versus sparse solver cost on growing tree sizes.
- Track nonlinear iteration counts and convergence behavior, not only runtime.

### Parallelization Strategy

The current tree solver is serial-only. Parallelization should be considered only after the serial structured-block path is correct and profiled.

Potential directions:

- Keep a serial tree solve on one rank initially for small/medium trees.
- Partition subtrees and perform local bottom-up condensation per rank.
- Communicate condensed subtree relations at partition boundaries.
- Recover top-down corrections with parent-to-child communication.
- Keep `NewtonSparse` or `Nox` as fallback for unsupported distributed layouts.

### Default Solver Policy

Recommended policy:

- Keep `Nox` as default until `NewtonTree` is stable across representative models and inputs.
- Use `NewtonSparse` as an intermediate debugging/reference path.
- Make `NewtonTree` the default only after broad validation and performance evidence.
- Remove NOX only if the replacement path fully covers required functionality and a fallback strategy exists.
