# Optimization Step 16 - Specialize Common Top-Down Recovery Shapes

## Scope

This step reduces serial `TreeNewtonLinearSolver` top-down recovery overhead for common grouped `2x2` element shapes.

Only the serial tree solver implementation was changed. The distributed tree solver was intentionally left unchanged.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Added Direct 2x2 Group Top-Down Recovery

The serial solver now has a specialized grouped top-down path for `block_size == 2` and child counts `0`, `1`, or `2`.

The new path loops once over the group and directly:

- Reads the element inlet-pressure correction.
- Writes the inlet pressure to `delta` through the Step 15 local-index path.
- Computes the two unknown correction values.
- Writes both unknown corrections directly to `delta`.
- Computes child inlet-pressure corrections directly when the element has children.

### Avoided Temporary Top-Down Buffers For Common Shapes

The specialized path avoids the generic grouped top-down staging arrays:

- `top_down_inlet_pressure_`
- `top_down_unknown_values_`
- `top_down_outlet_pressure_`
- `top_down_child_pressure_`

These buffers are still used by the generic fallback path for uncommon shapes.

### Kept Existing Dispatch Fallbacks

The existing scalar small-group path is unchanged.

The new specialized path is used only when:

```text
group_size > top_down_scalar_group_threshold
group.block_size == 2
group.child_count <= 2
```

All other grouped shapes still use `recover_top_down_group(...)`.

## Behavior Preserved

- The Newton equations are unchanged.
- The bottom-up condensation mathematics are unchanged.
- The top-down correction formulas are unchanged.
- The same inlet-pressure propagation is used for child elements.
- Correction values are still written through the Step 15 local-index delta path.
- The scalar small-tree path is unchanged.
- The generic grouped fallback path remains available.
- The distributed tree solver was not changed.

## Expected Performance Effect

This step targets `tree_top_down_s` by removing gather/compute/scatter staging work for common `2x2` groups.

The largest expected benefit is in larger serial trees where many same-shape `2x2` elements are recovered in grouped top-down traversal.

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

When benchmarking later, compare against Step 15 and Step 13 with focus on:

- `tree_top_down_s`
- `tree_solve_s`
- `linear_solve_s`
- `newton_total_s`
- nonlinear iteration count

## Conclusion

Step 16 keeps the same serial tree-solver behavior but adds a direct top-down grouped path for the common `2x2` shapes. It should be kept if release benchmarks show lower top-down recovery time without changing nonlinear convergence.
