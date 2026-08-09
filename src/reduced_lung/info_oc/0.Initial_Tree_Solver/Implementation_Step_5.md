# Implementation Step 5

Date: 2026-07-17

## Scope

Implemented a reduced-lung-specific linear solver abstraction for Newton correction systems and refactored the custom Newton solver to use it.

The baseline implementation of this abstraction still calls the existing sparse `Core::LinAlg::Solver`, so the numerical backend remains the same as before this step. No tree-based linear solver, tree metadata, input-file solver selector, residual equation, Jacobian equation, dof ordering, or row ordering was changed.

The `src/reduced_lung/src/1d_pipe_flow/` implementation was intentionally ignored.

## Files Changed

### `src/reduced_lung/src/4C_reduced_lung_linear_solver.hpp`

Added the reduced-lung Newton linear solver abstraction.

New types:

- `NewtonLinearSystemMetadata`
- `NewtonLinearSolver`
- `SparseNewtonLinearSolverContext`
- `SparseNewtonLinearSolver`

`NewtonLinearSolver` exposes one operation:

```text
solve(jacobian, residual, x, metadata, delta)
```

The sign convention is documented in the interface:

```text
jacobian * delta = -residual
x_new = x_old + delta
```

`NewtonLinearSystemMetadata` currently carries:

- current physical time,
- time-step size,
- nonlinear iteration index.

This metadata is intentionally small for now. It is the future extension point for tree topology, block layout, junction data, boundary-condition metadata, and diagnostics.

### `src/reduced_lung/src/4C_reduced_lung_linear_solver.cpp`

Implemented `SparseNewtonLinearSolver`.

Behavior:

- owns a `Core::LinAlg::Solver`,
- owns a right-hand-side vector on the Newton correction map,
- forms `rhs = -residual`,
- zeros `delta`,
- calls `Core::LinAlg::Solver::solve(jacobian, delta, rhs, solver_params)`,
- uses `refactor = true`,
- uses `reset = metadata.nonlinear_iteration == 0`,
- passes through an optional configured `Projector`,
- throws on nonzero linear solver status.

This preserves the current sparse-solver behavior while moving it behind a reduced-lung interface.

### `src/reduced_lung/src/4C_reduced_lung_newton_solver.hpp`

Refactored `NewtonSolverContext`.

Before this step, `NewtonSolver` directly received enough data to construct a `Core::LinAlg::Solver`.

After this step, `NewtonSolverContext` receives:

```text
std::shared_ptr<NewtonLinearSolver> linear_solver
```

The Newton solver now depends only on the reduced-lung linear-solver abstraction.

### `src/reduced_lung/src/4C_reduced_lung_newton_solver.cpp`

Refactored the Newton correction solve.

Before this step, `NewtonSolver::solve_linear_correction()` directly:

- formed the sparse RHS,
- configured `Core::LinAlg::SolverParams`,
- called `Core::LinAlg::Solver`.

After this step, it builds metadata and delegates to the interface:

```text
linear_solver_->solve(jacobian_, residual_, x_solution_, metadata, delta_)
```

The nonlinear loop remains unchanged:

- synchronize trial solution,
- assemble residual,
- check convergence,
- assemble Jacobian,
- solve for `delta`,
- update `x += delta`.

### `src/reduced_lung/tests/4C_reduced_lung_newton_solver_test.np2.cpp`

Updated the direct custom Newton solver test to construct a `SparseNewtonLinearSolver` and pass it to `NewtonSolverContext`.

The test behavior was otherwise unchanged.

### `src/reduced_lung/tests/4C_reduced_lung_newton_vs_nox_test.np2.cpp`

Updated the Newton-vs-NOX comparison fixture to construct a `SparseNewtonLinearSolver` and pass it to `NewtonSolverContext`.

The comparison coverage was otherwise unchanged.

## Preserved Behavior

- The existing `NoxSolver` implementation was not changed.
- `ReducedLungSimulation::solve_timestep()` still calls NOX in the production runtime path.
- No YAML input behavior was changed.
- No residual or Jacobian assembly code was changed.
- No dof or row ordering was changed.
- The custom Newton solver still uses the existing sparse linear solver backend through `SparseNewtonLinearSolver`.

## Purpose For Future Work

The custom Newton solver now has a clean insertion point for a future tree-based linear solver.

Current structure:

```text
NewtonSolver
  -> NewtonLinearSolver interface
      -> SparseNewtonLinearSolver
          -> Core::LinAlg::Solver / UMFPACK
```

Future tree-solver structure:

```text
NewtonSolver
  -> NewtonLinearSolver interface
      -> TreeNewtonLinearSolver
          -> tree metadata and bottom-up/top-down elimination
```

The nonlinear loop should not need to change when the tree linear solver is added.

## Verification Performed

Built the two-process reduced-lung unit-test executable:

- `cmake --build build/debug --target unittests_reduced_lung.np2 --parallel 4`

Ran all reduced-lung unit-test targets:

- `ctest -R "^unittests_reduced_lung(\.np2)?$" --output-on-failure`

Result: passed, 2 of 2 reduced-lung unit-test targets.

Ran representative non-1d reduced-lung input tests to confirm the unchanged NOX runtime path:

- `reduced_lung_serial_airways_flow.4C.yaml-p1`
- `reduced_lung_serial_airways_flow.4C.yaml-p2`
- `reduced_lung_aw_bifurcation_flow.4C.yaml-p1`
- `reduced_lung_aw_bifurcation_flow.4C.yaml-p2`
- `reduced_lung_3_aw_2_tu.4C.yaml-p1`
- `reduced_lung_3_aw_2_tu.4C.yaml-p3`
- `reduced_lung_terminal_unit.4C.yaml-p1`

Result: passed, 7 of 7 selected reduced-lung input tests. CTest also ran the `test_cleanup` fixture, which passed.

Attempted formatting:

- `clang-format` was not available on the current PATH, so no formatting command was applied.

## Current Limitations

- Only the sparse wrapper implementation exists.
- No tree-based linear solver has been added yet.
- No tree metadata has been added yet.
- The custom Newton solver is still not selectable from YAML input files.
