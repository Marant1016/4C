# Implementation Step 4

Date: 2026-07-17

## Scope

Implemented direct comparison tests between the existing reduced-lung NOX solver path and the new custom reduced-lung Newton solver path.

The comparison still uses the existing sparse linear solver backend, typically UMFPACK. No tree-based linear solver, tree metadata, input-file solver selector, residual equation, Jacobian equation, dof ordering, or row ordering was changed.

The `src/reduced_lung/src/1d_pipe_flow/` implementation was intentionally ignored.

## Files Changed

### `src/reduced_lung/tests/4C_reduced_lung_newton_vs_nox_test.np2.cpp`

Added a new two-process GoogleTest file dedicated to NOX-vs-custom-Newton comparisons.

This is a separate test file instead of an edit to `4C_reduced_lung_newton_solver_test.np2.cpp` because the existing Newton test is a direct analytical smoke test, while Phase 4 needs solver-to-solver comparison coverage across multiple model topologies.

The new test creates two independent but identical reduced-lung fixtures for each case:

- one fixture solved with `NoxSolver`,
- one fixture solved with `NewtonSolver`.

Each fixture builds the same reduced-lung objects used by the production setup path:

- `ReducedLungParameters`,
- local `Core::FE::Discretization`,
- airway and terminal-unit model containers,
- boundary-condition and junction containers,
- dof, row, and locally relevant maps,
- dof, locally relevant dof, and solution vectors,
- sparse Jacobian matrix,
- `ReducedLungAssemblyPipeline`,
- `Core::LinAlg::Solver` configuration using UMFPACK.

## Comparison Cases Added

### Single Terminal Unit With Ogden Elasticity

- Reuses the same analytical pressure function used by the existing NOX terminal-unit test.
- Compares NOX and custom Newton solution vectors after every time step.
- Compares owned dof vectors after every time step.
- Checks both residual norms are converged.
- Compares terminal-unit volume histories after end-of-timestep updates.

### Single Terminal Unit With Linear Elasticity

- Covers the second terminal-unit elasticity model.
- Uses Kelvin-Voigt rheology with linear elasticity.
- Compares solution vectors, owned dofs, residual convergence, and terminal-unit volume histories.

### Serial Rigid Airways

- Covers a serial airway topology with three rigid linear-resistive airway elements.
- Uses pressure boundary conditions at the inlet and outlet.
- Compares solution vectors and owned dofs.
- Checks both residual norms are converged.
- Checks connection flow balances for both solver results.

### Rigid-Airway Bifurcation

- Covers a one-parent/two-child bifurcation topology with rigid linear-resistive airway elements.
- Uses pressure boundary conditions at the inlet and both outlets.
- Compares solution vectors and owned dofs.
- Checks both residual norms are converged.
- Checks bifurcation flow balance for both solver results.

## Important Test Detail

The test compares the unique solution and owned dof vectors. It does not compare raw locally relevant/ghosted vectors across independently built fixtures because those column maps can contain non-unique global ids and can be distributed differently across MPI ranks.

The locally relevant vectors are still used inside each fixture for residual assembly, state updates, and flow-balance checks.

## Preserved Behavior

- The existing `NoxSolver` implementation was not changed.
- The existing `NewtonSolver` implementation was not changed for this phase.
- `ReducedLungSimulation::solve_timestep()` still calls NOX in the production runtime path.
- No YAML input behavior was changed.
- No residual or Jacobian assembly code was changed.
- No dof or row ordering was changed.

## Verification Performed

Built the two-process reduced-lung unit-test executable:

- `cmake --build build/debug --target unittests_reduced_lung.np2 --parallel 4`

Ran the two-process reduced-lung unit tests containing the new comparison tests:

- `ctest -R "^unittests_reduced_lung\.np2$" --output-on-failure`

Result: passed.

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

## Current Limitations

- The custom Newton solver is still not selectable from YAML input files.
- The comparison tests are C++ unit tests, not main-executable input-file tests.
- Exact nonlinear iteration histories are not required to match because NOX and the custom Newton loop perform convergence checks through different control code.
- The tests check converged residual norms and compare final solution values, but they do not expose or compare a full per-iteration increment history.
