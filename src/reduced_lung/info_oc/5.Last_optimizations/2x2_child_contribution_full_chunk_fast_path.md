# 2x2 Child Contribution Full-Chunk Fast Path

## Implementation Rules

When asked to implement a step:
1. Implement only that step.
2. Verify it with relevant tests.
3. Write a short summary to `src/reduced_lung/info_oc/5.Last_optimizations/`.

## Goal

Reduce overhead in `TreeNewtonLinearSolver::solve()` for the structured direct `2x2` bottom-up path.

Callgrind shows that `add_2x2_child_contribution_structured_chunk(...)` is dominated more by SIMD gather/scatter validity handling and indirect indexing than by the arithmetic itself.

The lowest-risk optimization is to split full SIMD chunks from the final partial chunk:

```text
full chunks:   gather(...) + scatter(...)
tail chunk:    existing gather_or(...) + scatter_valid(...)
```

This keeps current behavior for partial tails and avoids per-lane `grouped_index < valid_end` checks for full chunks.

## Target File

```text
src/reduced_lung/src/solver/tree/4C_reduced_lung_tree_linear_solver.cpp
```

## Target Functions And Lines

Primary target:
## Implementation Rules

When asked to implement a step:
1. Implement only that step.
2. Verify it with relevant tests.
3. Write a short summary to `src/reduced_lung/info_oc/5.Last_optimizations/`.
```text
add_2x2_child_contribution_structured_chunk(...)
```

Current location:

```text
approximately lines 1990-2083
```

Related helper functions:

```text
tree_solver_simd::gather(...)
tree_solver_simd::gather_or(...)
tree_solver_simd::scatter(...)
tree_solver_simd::scatter_valid(...)
tree_solver_simd::full_chunk_end(...)
tree_solver_simd::padded_chunk_end(...)
```

Current location:

```text
approximately lines 58-116
```

Call sites:

```text
assemble_2x2_one_child_group(...)
assemble_2x2_two_child_group(...)
```

Current location:

```text
approximately lines 2351-2402
```

Also update:

```text
subtract_2x2_rhs_shift_structured_chunk(...)
```

Current location:

```text
approximately lines 2086-2096
```

## Current Access Pattern

### Contiguous Accesses

These are contiguous within a grouped `2x2` batch:

```text
batch_2x2_a10_[lane_index]
batch_2x2_a11_[lane_index]
batch_2x2_rhs_constant1_[lane_index]
```

`lane_index` is computed as:

```text
grouped_index - group_begin
```

These are safe candidates for full-chunk `scatter(...)` instead of `scatter_valid(...)`.

### Indirect Accesses

These remain indirect and should not be converted to contiguous SIMD loads without a separate proof:

```text
grouped_child_begin_[grouped_index] + child_slot
pressure_row_[child_index]
child_element_index_[child_index]
subtree_relation_G_[child_element_index]
subtree_relation_h_[child_element_index]
parent_outlet_pressure_unknown_index_[child_index]
child_pressure_parent_coefficient_values_[child_index]
child_pressure_child_coefficient_values_[child_index]
child_flow_coefficient_values_[child_index]
```

Because of these indirect chains, do not assume contiguous SIMD loads are possible for this function.

## Step 1: Preserve Existing Tail-Safe Logic

Before changing behavior, keep the current implementation as the tail-safe path.

Suggested structure:

```text
add_2x2_child_contribution_structured_chunk(..., bool full_chunk)
```

or two small local lambdas inside the existing function:

```text
gather_child_value(...)
gather_child_value_or(...)
scatter_child_value(...)
scatter_child_value_valid(...)
```

Prefer the smallest readable change. Avoid adding a new named helper unless the lambda version becomes hard to follow.

## Step 2: Add Full-Chunk Variant Inside The Existing Function

For full chunks, replace:

```text
tree_solver_simd::gather_or(chunk_begin, valid_end, pad, load)
tree_solver_simd::scatter_valid(values, chunk_begin, valid_end, store)
```

with:

```text
tree_solver_simd::gather(chunk_begin, load)
tree_solver_simd::scatter(values, chunk_begin, store)
```

Only do this when the caller guarantees:

```text
chunk_begin + tree_solver_simd::width() <= group.end
```

Keep the current `gather_or(...)` and `scatter_valid(...)` implementation for the final partial chunk.

## Step 3: Split The 2x2 One-Child Loop

In `assemble_2x2_one_child_group(...)`, replace the current padded loop with two loops:

```text
const int full_end = tree_solver_simd::full_chunk_end(group.begin, group.end);

for (; grouped_index < full_end; grouped_index += tree_solver_simd::width())
{
  assemble_2x2_equation_rows_structured_chunk_full(...);
  tree_solver_simd::Double rhs_shift(0.0);
  add_2x2_child_contribution_structured_chunk_full(...);
  subtract_2x2_rhs_shift_structured_chunk_full(...);
}

if (grouped_index < group.end)
{
  assemble_2x2_equation_rows_structured_chunk_tail(...);
  tree_solver_simd::Double rhs_shift(0.0);
  add_2x2_child_contribution_structured_chunk_tail(...);
  subtract_2x2_rhs_shift_structured_chunk_tail(...);
}
```

The function names above are descriptive placeholders. The actual implementation can use a `bool full_chunk` parameter if that is simpler.

## Step 4: Split The 2x2 Two-Child Loop

Apply the same full/tail split in `assemble_2x2_two_child_group(...)`.

Full chunk:

```text
assemble equation rows
rhs_shift = 0
add child slot 0 with full path
add child slot 1 with full path
subtract rhs shift with full path
```

Tail chunk:

```text
assemble equation rows
rhs_shift = 0
add child slot 0 with current valid path
add child slot 1 with current valid path
subtract rhs shift with current valid path
```

## Step 5: Consider Equation Row Assembly Separately

`assemble_2x2_equation_rows_structured_chunk(...)` also uses `gather_or(...)` and `scatter_valid(...)`.

Do not change it in the first patch unless the `add_2x2_child_contribution_structured_chunk(...)` change is clean and measurable.

If the first patch helps, repeat the same full/tail split for equation-row assembly as a second independent patch.

## Step 6: Do Not Add Symbolic Caches In The First Patch

Repeated chains such as:

```text
grouped_child_begin_[grouped_index] + child_slot
child_element_index_[child_index]
pressure_row_[child_index]
parent_outlet_pressure_unknown_index_[child_index]
```

can be precomputed in the symbolic plan, but that is a larger change.

Defer it until after the full-chunk fast path is profiled.

If still needed, add narrowly scoped grouped child caches such as:

```text
grouped_child0_interface_index_
grouped_child1_interface_index_
grouped_child0_element_index_
grouped_child1_element_index_
grouped_child0_pressure_row_
grouped_child1_pressure_row_
grouped_child0_parent_outlet_pressure_unknown_index_
grouped_child1_parent_outlet_pressure_unknown_index_
```

Only add caches that remove measured hot indirections.

## Step 7: Preserve Correctness And Fallback Behavior

Do not change:

```text
require_structured_coefficients(...)
rhs_shift accumulation
child_pressure_slope_ writes
child_pressure_intercept_ writes
batch_2x2_a10_ / batch_2x2_a11_ update logic
tail handling
fallback behavior in solve_2x2_batch(...)
```

Full chunks must produce exactly the same values as the current padded valid-path implementation.

## Step 8: Verification Commands

Build:

```bash
cmake --build /scratch/Rodriguez/workspace/4C/4C/build/release \
  --target 4C benchmarktests_reduced_lung unittests_reduced_lung \
  --parallel 4
```

Run unit tests:

```bash
ctest -R "^unittests_reduced_lung$" --output-on-failure
```

Run nonlinear gen16 profile:

```bash
FOUR_C_REDUCED_LUNG_TREE_PROFILE=1 ./build/release/4C \
  ../files/gen16_inputs/reduced_lung_lung_tree_gen16_nonlinear_500steps_newton_tree.4C.yaml \
  ../output/gen16_nonlinear_newton_tree_profile
```

Run linear gen16 regression profile:

```bash
FOUR_C_REDUCED_LUNG_TREE_PROFILE=1 ./build/release/4C \
  ../files/reduced_lung_lung_tree_gen16_500steps_newton_tree.4C.yaml \
  ../output/gen16_newton_tree_profile
```

## Step 9: Metrics To Compare

Primary timing counters:

```text
tree_solve_s
tree_bottom_up_s
tree_dense_s
```

Correctness/direct-path invariants:

```text
tree_assembly_solver_update_s: 0
tree_lookup_s: 0
tree_lookups: 0
tree_dense_fallbacks: 0
tree_unsupported_fallbacks: 0
```

SIMD counters:

```text
tree_simd_groups
tree_simd_lanes
tree_scalar_tail_lanes
```

Expected result:

```text
tree_simd_groups and tree_simd_lanes should stay equivalent.
tree_dense_fallbacks should not increase.
tree_unsupported_fallbacks should not increase.
tree_bottom_up_s should be neutral or lower.
```

## Keep Or Revert Rule

Keep the change only if all of these hold:

```text
1. Unit tests pass.
2. Direct-vs-generic coefficient comparisons pass.
3. Sparse and direct tree corrections still match.
4. Direct-path invariants remain zero for lookup/update counters.
5. tree_bottom_up_s or Callgrind instruction count improves on nonlinear gen16.
6. Code remains readable and the tail-safe path is still obvious.
```

If timings are neutral and code readability worsens, do not keep the change.
