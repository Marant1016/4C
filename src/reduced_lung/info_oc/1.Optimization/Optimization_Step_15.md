# Optimization Step 15 - Write Delta Through Precomputed Local Indices

## Scope

This step reduces serial `TreeNewtonLinearSolver` top-down recovery overhead by replacing per-dof global correction writes with direct local vector writes.

Only the serial tree solver was changed. The distributed `DistributedTreeNewtonLinearSolver` still uses its existing correction-write path and was intentionally left unchanged.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Added Local Correction Write Indices

The serial solver now stores local correction indices for the dofs written during top-down recovery:

```cpp
std::vector<int> inlet_pressure_correction_local_dof_ids_;
std::vector<int> unknown_correction_local_dof_ids_;
bool correction_local_dof_ids_initialized_ = false;
```

These arrays are initialized once from `delta.get_map().lid(...)` during the first serial solve call. Each global correction dof is validated as locally available before the cache is marked initialized.

### Removed Serial Global Delta Writes

The old serial helper used:

```cpp
delta.replace_global_value(global_dof_id, value);
```

The serial top-down path now writes directly to the local vector storage:

```cpp
delta_values[local_dof_id] = value;
```

This avoids repeated Epetra global-id replacement work for every recovered correction dof.

### Updated Top-Down Recovery

Both serial top-down paths now use cached local indices:

- Scalar recovery through `recover_top_down_element(...)`
- Grouped recovery through `recover_top_down_group(...)`

The computed correction values are unchanged; only the vector write mechanism changed.

## Behavior Preserved

- The Newton equations are unchanged.
- The bottom-up tree condensation mathematics are unchanged.
- The top-down recovery formulas are unchanged.
- The same correction values are written to the same global dofs.
- Serial solves still assert that all correction dofs are locally available.
- Sparse-Jacobian and structured-tree coefficient behavior are unchanged.
- The distributed tree solver was not changed.

## Expected Performance Effect

This step targets `tree_top_down_s` by removing repeated global-id replacement calls in the serial correction scatter.

The largest expected benefit is in larger serial trees where top-down recovery writes many correction values.

## Verification

Verification passed:

```text
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
git diff --check
```

The `ctest` commands were run from `build/debug`.

Results:

- `unittests_reduced_lung`: passed, `1/1` tests.
- Focused `newton_tree` input tests: passed, `7/7` tests.
- `git diff --check`: passed.

## Benchmark Status

Release benchmarks were not run for this step.

When benchmarking later, compare against Step 14 and Step 13 with focus on:

- `tree_top_down_s`
- `tree_solve_s`
- `linear_solve_s`
- `newton_total_s`
- nonlinear iteration count

## Conclusion

Step 15 keeps the same serial tree-solver behavior but replaces repeated global correction writes with direct local writes. It should be kept if release benchmarks show lower top-down recovery time without changing nonlinear convergence.
