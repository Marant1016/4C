# Optimization Step 12 - Specialize Common Element Shapes

## Scope

This step specializes common serial tree bottom-up assembly shapes.

Only the serial `TreeNewtonLinearSolver` bottom-up assembly path was changed. The distributed `DistributedTreeNewtonLinearSolver` implementation was intentionally left unchanged.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Added Specialized `2x2` Assembly Helpers

The serial solver now has specialized bottom-up assembly paths for common `2x2` element groups:

```text
block_size == 2, child_count == 0
block_size == 2, child_count == 1
block_size == 2, child_count == 2
```

These are handled by small fixed-shape helpers that assemble the two local equation rows directly and avoid the generic equation loop.

### Specialized Child Contributions

The one-child and two-child `2x2` paths compute child pressure and flow contributions directly for the fixed child slots.

The formula is unchanged:

```cpp
child_pressure_slope = -pressure_parent_coeff / pressure_child_coeff;
child_pressure_intercept = pressure_rhs / pressure_child_coeff;
child_flow_slope = subtree_relation_G[child] * child_pressure_slope;
child_flow_intercept = subtree_relation_G[child] * child_pressure_intercept + subtree_relation_h[child];
```

The resulting child-flow contribution is still added to the parent outlet-pressure matrix entry, and the RHS shift is still subtracted from the last local equation row.

### Kept Generic Fallback

Bottom-up assembly now dispatches common `2x2` groups to the specialized paths. All other shapes still use the existing generic `assemble_group(...)` implementation.

The existing dense solve dispatch is unchanged: `2x2` groups still use `solve_2x2_batch(...)`, `3x3` groups still use `solve_3x3_batch(...)`, and uncommon sizes still use the scalar dense solver.

## Behavior Preserved

- The bottom-up condensation mathematics is unchanged.
- The top-down recovery path from Step 11 is unchanged.
- Direct coefficient storage from Step 10 is unchanged.
- The generic assembly fallback remains available for uncommon shapes.
- Required coefficient assertions are unchanged.
- Sparse-Jacobian fallback behavior is unchanged.
- The distributed tree solver was not changed.

## Intended Benefit

The generic bottom-up assembly loop handles all shapes and child counts, which adds loop and branch overhead to the common `2x2` cases. This step removes that overhead for the common leaf, one-child, and two-child element groups while keeping the existing grouped dense solver path.

The expected benchmark signal is lower `tree_bottom_up_s`, especially in cases where `tree_max_block=2` and most groups fall into the specialized shapes.

## Verification

Verification passed:

```text
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

The reduced-lung unit test and focused `newton_tree` input tests passed.

No release benchmark rerun was performed in this step.

## Follow-Up

The next release benchmark run should compare `tree_bottom_up_s`, `tree_solve_s`, and full-solve `newton_total_s` against Steps 10-11. If this specialization does not produce a measurable benefit, the helper duplication should be reconsidered before adding more shape-specific kernels.
