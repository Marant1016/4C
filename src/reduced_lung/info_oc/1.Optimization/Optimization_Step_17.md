# Optimization Step 17 - Remove Redundant Batch Dense-Solver Staging

## Scope

This step reduces serial `TreeNewtonLinearSolver` bottom-up dense-solve overhead for grouped `2x2` local blocks.

Only the serial tree solver implementation and its serial solver workspace declarations were changed. The distributed tree solver was intentionally left unchanged.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Removed 2x2 Staging Arrays

The serial solver no longer stores the temporary structure-of-arrays buffers used only by the old `2x2` batch solve path:

```cpp
batch_2x2_a00_
batch_2x2_a01_
batch_2x2_a10_
batch_2x2_a11_
batch_2x2_rhs_constant0_
batch_2x2_rhs_constant1_
batch_2x2_rhs_inlet_pressure0_
batch_2x2_rhs_inlet_pressure1_
batch_2x2_intercept0_
batch_2x2_intercept1_
batch_2x2_slope0_
batch_2x2_slope1_
```

Only `batch_2x2_fallback_lanes_` remains because it is still needed to preserve the dense fallback for singular or near-singular `2x2` blocks.

### Solved 2x2 Blocks Directly From Existing Workspaces

`solve_2x2_batch(...)` now reads directly from:

- `workspace_matrix_`
- `workspace_rhs_constant_`
- `workspace_rhs_inlet_pressure_`

It writes directly into:

- `workspace_intercept_`
- `workspace_slope_`

This removes the previous copy from generic workspace arrays into the temporary `batch_2x2_*` arrays before solving.

### Preserved Dense Fallback

If a determinant is not finite or is below `pivot_tolerance`, the solver still records the lane and calls the existing generic dense solver fallback:

```cpp
solve_dense_system(..., 2, pivot_tolerance, ...)
```

## Behavior Preserved

- The local `2x2` equation math is unchanged.
- The intercept and slope formulas are unchanged.
- Near-singular fallback behavior is unchanged.
- Bottom-up tree condensation results are unchanged.
- `3x3` batch solving is unchanged.
- The distributed tree solver was not changed.

## Expected Performance Effect

This step targets `tree_bottom_up_s` and `tree_dense_s` by removing memory traffic and staging overhead in the common grouped `2x2` dense-solve path.

The largest expected benefit is in larger structured serial trees with many grouped `2x2` element solves.

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

When benchmarking later, compare against Step 16 and Step 13 with focus on:

- `tree_bottom_up_s`
- `tree_dense_s`
- `tree_solve_s`
- `linear_solve_s`
- `newton_total_s`
- nonlinear iteration count

## Conclusion

Step 17 keeps the same solver behavior but removes redundant staging in the common `2x2` batch dense-solve path. It should be kept if release benchmarks show lower bottom-up or dense-solve time without changing nonlinear convergence.
