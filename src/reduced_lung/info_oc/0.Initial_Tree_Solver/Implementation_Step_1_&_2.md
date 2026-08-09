# Implementation Step 1 & 2

Date: 2026-07-17

## Scope

Implemented preparation work for the first two tree-based solver plan phases:

1. Preserve baseline NOX behavior.
2. Generalize the reduced-lung assembly pipeline so it can be reused by a future custom Newton solver.

The `src/reduced_lung/src/1d_pipe_flow/` implementation was intentionally ignored.

## Files Changed

### `src/reduced_lung/src/4C_reduced_lung_helpers.hpp`

- Introduced `ReducedLungAssemblyPipeline` as the neutral assembly callback registry name.
- Kept the same callback groups:
  - `residual_assemblers`
  - `jacobian_assemblers`
  - `state_updaters`
- Added `using NoxAssemblyPipeline = ReducedLungAssemblyPipeline;` as a compatibility alias for existing NOX-oriented code and tests.
- Added `create_default_reduced_lung_assembly_pipeline(...)` as the neutral factory declaration.
- Kept `create_default_nox_assembly_pipeline(...)` as a compatibility factory declaration.
- Updated `NoxSolverContext::assembly_pipeline` and `NoxSolver::assembly_pipeline_` to use `ReducedLungAssemblyPipeline`.

### `src/reduced_lung/src/4C_reduced_lung_helpers.cpp`

- Renamed the default assembly pipeline implementation to `create_default_reduced_lung_assembly_pipeline(...)`.
- Preserved the existing callback order exactly:
  - Residual: Airways, Terminal units, Junctions, Boundary conditions.
  - Jacobian: Airways, Terminal units, Junctions, Boundary conditions.
  - State updates: Airways, Terminal units.
- Added `create_default_nox_assembly_pipeline(...)` as a wrapper that calls `create_default_reduced_lung_assembly_pipeline(...)`.
- Did not change `NoxSolver::solve`, `NoxSolver::residual`, `NoxSolver::jacobian`, or `NoxSolver::sync_state_from_x` behavior.

### `src/reduced_lung/src/4C_reduced_lung_main.cpp`

- Updated `ReducedLungSimulation::build_linear_system_and_solver()` to create the assembly pipeline through `create_default_reduced_lung_assembly_pipeline(...)`.
- Updated the stored pipeline member type to `ReducedLungAssemblyPipeline`.
- Left NOX as the active runtime nonlinear solver.
- Left `ReducedLungSimulation::solve_timestep()` calling `nox_solver_->solve(current_time_)`.

## Preserved Behavior

- No residual equations were changed.
- No Jacobian equations were changed.
- No dof ordering was changed.
- No row ordering was changed.
- No input-file behavior or solver-selection behavior was changed.
- The current NOX path remains the only runtime path in `ReducedLungSimulation`.

## Purpose For Future Work

The reduced-lung assembly callbacks are now named independently of NOX. A future custom Newton solver can reuse `ReducedLungAssemblyPipeline` while the existing NOX solver continues to use the same callback registry through the compatibility alias and wrapper.

## Verification Performed

Built targeted debug unit-test executables:

- `cmake --build build/debug --target unittests_reduced_lung --parallel 4`
- `cmake --build build/debug --target unittests_reduced_lung.np2 --parallel 4`

Ran targeted reduced-lung unit tests:

- `ctest -R "^unittests_reduced_lung(\.np2)?$" --output-on-failure`

Result: passed, 2 of 2 tests.

Ran representative non-1d reduced-lung input tests:

- `reduced_lung_serial_airways_flow.4C.yaml-p1`
- `reduced_lung_serial_airways_flow.4C.yaml-p2`
- `reduced_lung_aw_bifurcation_flow.4C.yaml-p1`
- `reduced_lung_aw_bifurcation_flow.4C.yaml-p2`
- `reduced_lung_3_aw_2_tu.4C.yaml-p1`
- `reduced_lung_3_aw_2_tu.4C.yaml-p3`
- `reduced_lung_terminal_unit.4C.yaml-p1`

Result: passed, 7 of 7 selected reduced-lung input tests. CTest also ran the `test_cleanup` fixture, which passed.
