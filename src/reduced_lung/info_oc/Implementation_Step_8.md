# Implementation Step 8

Date: 2026-07-18

## Scope

Added runtime solver-workflow selection for the reduced-lung model.

The implementation keeps NOX as the default solver workflow and adds selectable custom Newton workflows using either the sparse linear solver wrapper or the tree linear solver.

NOX was intentionally not removed. It remains the production default and the reference/fallback path while the custom Newton and tree solver paths continue to be validated.

The `src/reduced_lung/src/1d_pipe_flow/` implementation was intentionally ignored.

## Files Changed

### `src/reduced_lung/src/4C_reduced_lung_input.hpp`

Added a nonlinear solver workflow enum:

```text
ReducedLungParameters::NonlinearSolverType
```

Supported values:

- `Nox`
- `NewtonSparse`
- `NewtonTree`

Added this field to `ReducedLungParameters::Dynamics`:

```text
nonlinear_solver
```

The default is:

```text
Nox
```

This preserves existing input-file behavior when the new field is omitted.

### `src/reduced_lung/src/4C_reduced_lung_input.cpp`

Registered the new input parameter in the `reduced_dimensional_lung.dynamics` group:

```text
nonlinear_solver
```

Description:

```text
Nonlinear solver workflow: Nox, NewtonSparse, or NewtonTree.
```

Default:

```text
Nox
```

Example future YAML usage:

```yaml
reduced_dimensional_lung:
  dynamics:
    nonlinear_solver: NewtonTree
```

### `src/reduced_lung/src/4C_reduced_lung_main.cpp`

Changed `ReducedLungSimulation::build_linear_system_and_solver()` so it constructs one selected nonlinear solver workflow instead of always constructing NOX.

Current selection behavior:

```text
Nox
  -> NoxSolver

NewtonSparse
  -> NewtonSolver
  -> SparseNewtonLinearSolver
  -> Core::LinAlg::Solver

NewtonTree
  -> build_reduced_lung_tree_metadata(...)
  -> NewtonSolver
  -> TreeNewtonLinearSolver
```

Added helper setup methods inside `ReducedLungSimulation`:

- `build_nox_solver()`
- `build_newton_solver_with_sparse_linear_solver()`
- `build_newton_solver_with_tree_linear_solver()`
- `build_newton_solver()`

Changed `ReducedLungSimulation::solve_timestep()` so it calls the selected solver:

```text
NoxSolver::solve(current_time)
```

or:

```text
NewtonSolver::solve(current_time)
```

The existing end-of-time-step routines are unchanged:

- `TerminalUnits::end_of_timestep_routine(...)`
- `Airways::end_of_timestep_routine(...)`

Added storage for the custom Newton workflow:

- `newton_solver_`
- `newton_linear_solver_`
- `tree_metadata_`

The `NewtonTree` path currently checks for serial execution before building metadata because the first tree solver is serial-only.

### `src/reduced_lung/tests/4C_reduced_lung_input_pipeline_test.cpp`

Added a unit test confirming that the new nonlinear solver workflow defaults to `Nox`.

## Solver Workflow Meanings

### `Nox`

Current production behavior.

```text
NOX nonlinear loop
  -> existing residual assembly
  -> existing Jacobian assembly
  -> Core::LinAlg::Solver / UMFPACK
```

### `NewtonSparse`

Custom Newton loop using the existing sparse linear solver backend.

```text
NewtonSolver
  -> existing residual assembly
  -> existing Jacobian assembly
  -> SparseNewtonLinearSolver
  -> Core::LinAlg::Solver / UMFPACK
```

This path is useful to isolate custom Newton behavior from tree-solver behavior.

### `NewtonTree`

Custom Newton loop using the first tree linear solver implementation.

```text
NewtonSolver
  -> existing residual assembly
  -> existing Jacobian assembly
  -> TreeNewtonLinearSolver
  -> bottom-up condensation / top-down recovery
```

This path currently still assembles the sparse Jacobian and uses it as the coefficient source for the tree solver. A future optimization can replace global sparse Jacobian assembly with structured derivative blocks.

## NOX Retention Decision

NOX was kept intentionally.

Reasons:

- NOX is the current production reference path.
- The tree solver is still serial-only.
- The tree solver has only been validated on selected small systems so far.
- NOX is useful for regression tests and fallback while the new path matures.
- Removing `solver_nonlin_nox` now would be premature.

Recommended future progression:

1. Keep `Nox` as the default while validating `NewtonSparse` and `NewtonTree`.
2. Use `NewtonSparse` to isolate custom Newton-loop behavior.
3. Use `NewtonTree` as the target tree-solver workflow.
4. Only make a custom path the default after representative input tests are stable.
5. Only remove the NOX dependency after the replacement path supports required models, topology cases, and execution modes.

## Preserved Behavior

- Existing input files that omit `nonlinear_solver` still use `Nox`.
- The production default remains NOX.
- Existing residual assembly was not changed.
- Existing Jacobian assembly was not changed.
- Dof ordering was not changed.
- Row ordering was not changed.
- The reduced-lung CMake dependency on `solver_nonlin_nox` remains.

## Verification Performed

Built the serial reduced-lung unit-test executable:

```text
cmake --build build/debug --target unittests_reduced_lung --parallel 4
```

Result: passed.

Built the two-process reduced-lung unit-test executable:

```text
cmake --build build/debug --target unittests_reduced_lung.np2 --parallel 4
```

Result: passed.

Ran serial and `.np2` reduced-lung unit tests:

```text
ctest -R "^unittests_reduced_lung(\.np2)?$" --output-on-failure
```

Result: passed, 2 of 2 selected test targets.

Ran representative reduced-lung runtime input tests that continue to use the default NOX path:

```text
ctest -R "^(reduced_lung_terminal_unit\.4C\.yaml-p1|reduced_lung_serial_airways_flow\.4C\.yaml-p1|reduced_lung_aw_bifurcation_flow\.4C\.yaml-p1)$" --output-on-failure
```

Result: passed, 3 of 3 selected input tests. CTest also ran the `test_cleanup` fixture, which passed.

## Current Limitations

- `NewtonTree` remains serial-only.
- `NewtonTree` still uses the assembled sparse Jacobian as its coefficient source.
- No existing YAML input files were changed to use `NewtonSparse` or `NewtonTree` by default.
- NOX remains a required reduced-lung module dependency because it is still selectable and still the default.
