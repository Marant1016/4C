# Implementation Step 3

Date: 2026-07-17

## Scope

Implemented the first custom reduced-lung Newton solver beside the existing NOX solver.

The implementation intentionally keeps the current NOX runtime path unchanged and continues to ignore `src/reduced_lung/src/1d_pipe_flow/`.

## Files Changed

### `src/reduced_lung/src/4C_reduced_lung_newton_solver.hpp`

- Added `ReducedLung::NewtonSolverContext`.
- Added `ReducedLung::NewtonSolver`.
- The context mirrors the objects needed by `NoxSolver`, but binds a concrete `Core::LinAlg::SparseMatrix&` because the custom Newton path assembles and solves sparse Newton systems directly.
- The solver stores references to the nonlinear solution vector, owned dof vector, locally relevant dof vector, Jacobian matrix, and a copy of `ReducedLungAssemblyPipeline`.
- The solver owns residual, right-hand-side, and correction vectors.

### `src/reduced_lung/src/4C_reduced_lung_newton_solver.cpp`

- Implemented the custom full-step Newton loop.
- Reuses `ReducedLungAssemblyPipeline` callbacks for:
  - state synchronization updates,
  - residual assembly,
  - Jacobian assembly.
- Uses the current sparse linear solver infrastructure through `Core::LinAlg::Solver`.
- Solves Newton correction systems with the sign convention:
  - `J * delta = -F`
  - `x = x + delta`
- Uses the same synchronization pattern as the NOX wrapper:
  - export nonlinear trial `x` to owned dofs,
  - export owned dofs to locally relevant dofs,
  - run registered state updater callbacks.
- Preserves the existing matrix fill behavior:
  - first Jacobian assembly inserts sparsity and calls `complete()` if needed,
  - later assemblies update the already-filled matrix without resetting it.
- Uses `Core::LinAlg::SolverParams` with:
  - `refactor = true`,
  - `reset = iteration == 0`,
  - optional `Projector` pass-through if configured on the linear solver.
- Throws on nonlinear non-convergence after `max_nonlinear_iterations` corrections.
- Throws on nonzero linear solver status.

### `src/reduced_lung/tests/4C_reduced_lung_newton_solver_test.np2.cpp`

- Added a focused two-process test for the custom Newton solver.
- The test mirrors the existing NOX single-terminal-unit Ogden setup.
- It verifies the analytical volume behavior `V(t) = 1 + t` over five time steps.
- It constructs `NewtonSolver` directly and does not change the production runtime path.

## Preserved Behavior

- The existing `NoxSolver` implementation was not changed.
- `ReducedLungSimulation::solve_timestep()` still calls `nox_solver_->solve(current_time_)`.
- No input-file option or solver-selection behavior was changed.
- No residual equations were changed.
- No Jacobian equations were changed.
- No dof ordering was changed.
- No row ordering was changed.

## Current Limitations

- The custom Newton solver is implemented and tested directly, but it is not yet selectable from the reduced-lung input file.
- The custom Newton solver still uses `Core::LinAlg::Solver`; no tree-based linear solver interface or tree metadata is introduced yet.
- The nonlinear strategy is full-step Newton only, matching the current NOX full-step line-search configuration at the mathematical update level but not reproducing all NOX status-test infrastructure.

## Verification Performed

Built the two-process reduced-lung unit-test executable:

- `cmake --build build/debug --target unittests_reduced_lung.np2 --parallel 4`

Ran targeted reduced-lung unit tests:

- `ctest -R "^unittests_reduced_lung\.np2$" --output-on-failure`
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
