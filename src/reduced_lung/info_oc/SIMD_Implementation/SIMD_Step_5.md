# SIMD Step 5

Date: 2026-08-03

## Goal

Add SIMD-aware `2x2` top-down recovery for the serial `TreeNewtonLinearSolver`.

This step targets the specialized top-down recovery path:

```text
recover_2x2_top_down_group(...)
```

It follows the previous bottom-up SIMD work from Steps 3 and 4 and targets `tree_top_down_s`.

## Files Changed

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
```

## Implementation

Refactored the existing scalar `2x2` top-down body into a lane helper:

```text
recover_2x2_top_down_lane(...)
```

Added a guarded SIMD chunk helper when `<experimental/simd>` is available:

```text
recover_2x2_top_down_chunk(...)
```

The group helper now follows this structure:

```text
SIMD full chunks
scalar tail
```

When SIMD support is unavailable, it uses the scalar lane helper for the whole group.

## SIMD Arithmetic

For each full SIMD chunk, the implementation gathers inlet pressures and `2x2` affine recovery data:

```text
inlet_pressure
slope0, slope1
intercept0, intercept1
```

Then it computes correction values with SIMD arithmetic:

```text
value0 = slope0 * inlet_pressure + intercept0
value1 = slope1 * inlet_pressure + intercept1
```

For groups with children, it also computes outlet pressure:

```text
outlet_pressure = outlet_slope * inlet_pressure + outlet_intercept
```

For each child slot, it computes child inlet pressure:

```text
child_pressure = child_pressure_slope * outlet_pressure + child_pressure_intercept
```

## Scatter Writes

The final writes remain scalar scatter operations because correction dof ids and child element ids are indirect.

The SIMD path still calls the same write helpers:

```text
set_delta_local_value(...)
set_inlet_pressure(...)
```

This preserves existing validation for missing inlet-pressure corrections and duplicate child inlet-pressure writes.

## Data Layout Used

The SIMD chunk path uses grouped caches from Step 2:

```text
grouped_unknown_begin_
grouped_child_begin_
grouped_element_indices_
```

It also reads the existing top-down recovery work arrays:

```text
workspace_slope_
workspace_intercept_
child_pressure_slope_
child_pressure_intercept_
inlet_pressure_by_element_
```

## Fallback Behavior

Scalar recovery remains active for:

- small trees using `use_scalar_tree_solve_`,
- groups with size `<= top_down_scalar_group_threshold`,
- scalar tails after full SIMD chunks,
- builds without `<experimental/simd>`,
- non-`2x2` groups,
- generic `recover_top_down_group(...)`.

Unlike the dense solve, no pivot fallback is needed because top-down recovery only evaluates already-computed affine relations.

## Behavior Preserved

- The small-tree scalar path is unchanged.
- The top-down grouped dispatch is unchanged.
- Existing correction write ownership is unchanged.
- `set_inlet_pressure(...)` validation remains active.
- `delta` values are still written through the same local-dof ids.
- Bottom-up assembly and dense solve paths are unchanged.
- `3x3` and generic top-down paths are unchanged.
- The distributed tree solver is unchanged.

## Expected Benefit

This step targets `tree_top_down_s` for large `2x2` groups.

The best expected benefit is in large groups with children, where both correction recovery and child inlet-pressure propagation are vectorized.

The speedup may be moderate because the final writes are still indirect scalar scatters.

## Verification

Verification target for this step:

```text
git diff --check
cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

Performance check target after a release build is available:

```text
./build/release/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/LinearSolve/BalancedAirways/StructuredTree
```

Record at least:

```text
tree_solve_s
tree_bottom_up_s
tree_top_down_s
tree_dense_s
```

## Follow-Up

After this step, the main explicit SIMD coverage for common serial structured `2x2` groups is:

```text
bottom-up assembly
2x2 dense solve
top-down recovery
```

The next SIMD candidate is the `3x3` dense batch solve if benchmark data shows enough `3x3` work to justify the added complexity.
