# SIMD Step New 7

Date: 2026-08-04

## Goal

Batch subtree relation writes for the serial `TreeNewtonLinearSolver` bottom-up path. Before this step, every bottom-up group wrote subtree relations through `write_subtree_relation_group(...)`, which looped over elements and read the interleaved `workspace_slope_` and `workspace_intercept_` arrays.

This step adds specialized writers that use the existing SoA batch solve outputs for supported production block sizes.

## Files Changed

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
```

## Specialized Writers

Added:

```text
write_2x2_subtree_relation_group(...)
write_3x3_subtree_relation_group(...)
```

The `2x2` writer reads from:

```text
batch_2x2_slope0_
batch_2x2_slope1_
batch_2x2_intercept0_
batch_2x2_intercept1_
```

The `3x3` writer reads from:

```text
batch_3x3_slope0_
batch_3x3_slope1_
batch_3x3_slope2_
batch_3x3_intercept0_
batch_3x3_intercept1_
batch_3x3_intercept2_
```

Both writers select the correct SoA lane value using `inlet_flow_unknown_index_` and write:

```text
subtree_relation_G_[element] = selected slope
subtree_relation_h_[element] = selected intercept
```

## SIMD Path

For SIMD-enabled builds, the new writers process groups in padded valid-lane chunks using:

```text
tree_solver_simd::padded_chunk_end(...)
tree_solver_simd::gather_or(...)
tree_solver_simd::scatter_valid(...)
```

Non-SIMD builds use scalar lane fallbacks inside the specialized writers.

## Dispatch Changes

The bottom-up group completion now dispatches subtree relation writes by block size:

```cpp
if (group.block_size == 2)
{
  write_2x2_subtree_relation_group(group);
}
else if (group.block_size == 3)
{
  write_3x3_subtree_relation_group(group);
}
else
{
  write_subtree_relation_group(group);
}
```

The generic writer remains available for scalar small-tree traversal and unsupported future block sizes.

## Behavior Preserved

- `workspace_intercept_` and `workspace_slope_` write-back remains in place because top-down recovery still reads those arrays.
- Sparse-Jacobian and non-direct batch paths are supported because they pack into SoA before solving.
- Dense fallback solves write their results back into the SoA outputs before subtree relation writing.
- Unsupported future block sizes still use the generic workspace-based writer.
- The distributed tree solver was not changed.

## Verification

Commands run after implementation:

```text
cmake --build build/debug --target unittests_reduced_lung --parallel 4
git diff --check
ctest -R "^unittests_reduced_lung$" --output-on-failure
```
