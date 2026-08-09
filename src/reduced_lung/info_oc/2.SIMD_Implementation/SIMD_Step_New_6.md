# SIMD Step New 6

Date: 2026-08-04

## Goal

Add a dedicated SIMD top-down recovery path for serial `TreeNewtonLinearSolver` `3x3` groups. Before this step, large `3x3` groups used the generic scalar `recover_top_down_group(...)` path during top-down correction recovery.

This step covers:

- `3x3 + 0 child`
- `3x3 + 1 child`
- `3x3 + 2 children`

The existing small top-down group scalar bypass remains unchanged.

## Files Changed

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
```

## New 3x3 Recovery Path

Added `recover_3x3_top_down_group(...)` beside the existing `recover_2x2_top_down_group(...)` implementation.

The new path validates:

- `group.block_size == 3`
- `0 <= group.child_count <= 2`

For SIMD-enabled builds, it processes groups in padded valid-lane chunks using `tree_solver_simd::padded_chunk_end(...)`, `gather_or(...)`, `scatter_valid(...)`, and `for_each_valid_lane(...)`.

## SIMD Work Per Chunk

Each chunk does the following:

- Checks that each valid element has an inlet-pressure stamp for the current solve.
- Writes each element inlet-pressure correction to `delta`.
- Gathers the three local slope values and three local intercept values from `workspace_slope_` and `workspace_intercept_`.
- Computes the three unknown corrections as SIMD values.
- Scatters the three unknown corrections to the corresponding local correction dofs.
- For child groups, gathers the outlet pressure slope/intercept using `outlet_pressure_unknown_index_` and computes outlet pressure.
- For each child slot, gathers child pressure slope/intercept, computes child inlet pressure, and calls `set_inlet_pressure(...)` for valid lanes.

## Scalar Fallback

The new recovery helper includes a scalar lane fallback for non-SIMD builds. It mirrors the SIMD work without changing the generic recovery path for unsupported shapes.

## Dispatch Changes

Top-down group dispatch now routes supported `3x3` groups through the dedicated recovery path:

```cpp
if (group.block_size == 2 && group.child_count <= 2)
{
  recover_2x2_top_down_group(group);
}
else if (group.block_size == 3 && group.child_count <= 2)
{
  recover_3x3_top_down_group(group);
}
else
{
  recover_top_down_group(group);
}
```

The existing `top_down_scalar_group_threshold` behavior is preserved, so tiny groups still use the scalar element recovery path.

## Behavior Preserved

- `workspace_slope_` and `workspace_intercept_` remain the source for top-down recovery.
- Inlet pressure stamp checks and duplicate-stamp protection are unchanged.
- Child inlet pressures are still written through `set_inlet_pressure(...)`.
- Generic scalar recovery remains available for unsupported future shapes.
- The distributed tree solver was not changed.

## Verification

Commands run after implementation:

```text
cmake --build build/debug --target unittests_reduced_lung --parallel 4
git diff --check
ctest -R "^unittests_reduced_lung$" --output-on-failure
```
