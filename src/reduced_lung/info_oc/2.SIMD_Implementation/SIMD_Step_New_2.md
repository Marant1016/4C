# SIMD Step New 2

Date: 2026-08-04

## Goal

Replace scalar tail handling in the existing serial `TreeNewtonLinearSolver` SIMD-enabled paths with padded SIMD chunks and valid-lane guarded stores.

This step only changes tail mechanics for paths that already had SIMD chunks. It does not add new 3x3 assembly kernels, does not change tree traversal order, and does not remove numerical dense-solve fallbacks.

## Files Changed

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
```

## SIMD Helpers Added

Added helper utilities in the internal `tree_solver_simd` namespace:

```cpp
padded_chunk_end(...)
gather_or(...)
scatter_valid(...)
for_each_valid_lane(...)
```

Purpose:

- `padded_chunk_end(...)` rounds a group range up to a full SIMD-width chunk.
- `gather_or(...)` loads real lanes and returns a safe pad value for invalid lanes.
- `scatter_valid(...)` stores only lanes whose grouped index is inside the real group range.
- `for_each_valid_lane(...)` executes callbacks only for valid lanes in a padded chunk.

## Batch Solve Tail Changes

Updated:

```text
solve_2x2_batch(...)
solve_3x3_batch(...)
```

Before this step:

- Full SIMD chunks were processed first.
- Remaining lanes were solved by scalar tail loops.
- `scalar_tail_lane_count` was incremented for these lanes.

After this step:

- The final partial group is processed as a padded SIMD chunk when SIMD is available.
- Invalid padded lanes use identity matrix entries and zero RHS values.
- Invalid padded lanes are never stored and are never added to fallback lists.
- `simd_lane_count` counts the real group size instead of only full-width chunks.

Padding values:

```text
2x2: identity matrix, zero RHS
3x3: identity matrix, zero RHS
```

This keeps padded lanes finite and pivot-safe without writing their results.

## Structured 2x2 Assembly Tail Changes

Updated the structured SIMD branches in:

```text
assemble_2x2_leaf_group(...)
assemble_2x2_one_child_group(...)
assemble_2x2_two_child_group(...)
```

Before this step:

- Full SIMD chunks used structured coefficient SIMD assembly.
- Remaining lanes dropped into scalar element assembly.

After this step:

- Structured groups loop to `padded_chunk_end(...)`.
- Equation-row assembly uses `gather_or(...)` and `scatter_valid(...)`.
- Child contribution assembly uses valid-lane guarded coefficient checks and stores.
- Scalar tails were removed from the structured SIMD branches.

Important detail:

`require_structured_coefficients(...)` now checks only valid lanes in a padded chunk. This prevents padded lanes from triggering coefficient assertions while preserving checks for real elements.

## 2x2 Top-Down Recovery Tail Changes

Updated:

```text
recover_2x2_top_down_group(...)
```

Before this step:

- Full SIMD chunks used `recover_2x2_top_down_chunk(...)`.
- Remaining lanes used scalar `recover_2x2_top_down_lane(...)`.

After this step:

- The SIMD branch loops to `padded_chunk_end(...)`.
- Inlet pressure checks and writes are guarded by `for_each_valid_lane(...)`.
- Correction writes and child inlet-pressure writes use `scatter_valid(...)`.
- Invalid padded lanes use zero arithmetic inputs and do not write anything.

The scalar lane lambda remains for non-SIMD builds and is marked `[[maybe_unused]]` to avoid warning-as-error failures when SIMD is available.

## Behavior Preserved

- Existing scalar paths remain for builds without `<experimental/simd>`.
- Existing small-tree scalar traversal remains unchanged.
- Existing generic `assemble_group(...)` remains unchanged.
- Existing generic `recover_top_down_group(...)` remains unchanged.
- Existing dense fallback logic remains unchanged for real lanes.
- Sparse-Jacobian validation behavior remains unchanged outside the affected SIMD-enabled paths.

## Expected Counter Impact

For SIMD-enabled large-tree runs, this step should reduce or eliminate scalar tail accounting from the updated paths:

```text
tree_scalar_tail_lanes == 0
```

Other scalar counters may still be nonzero because later steps are still required for generic 3x3 assembly/recovery, small-tree mode, unsupported block sizes, or numerical fallback cases.

## Verification

Commands run after this step:

```text
cmake --build build/debug --target unittests_reduced_lung --parallel 4
```

Recommended final checks:

```text
git diff --check
ctest -R "^unittests_reduced_lung$" --output-on-failure
```
