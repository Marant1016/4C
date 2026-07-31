# Optimization Step 8 - Vectorize Serial Top-Down Batches

## Scope

This step implements Phase 8 from the serial tree solver vectorization guide.

Only the serial `TreeNewtonLinearSolver` top-down recovery path was changed. The distributed `DistributedTreeNewtonLinearSolver` implementation was intentionally left unchanged.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Added Top-Down Batch Workspace

The serial solver now owns reusable lane buffers for grouped top-down recovery:

```cpp
std::vector<double> top_down_inlet_pressure_;
std::vector<double> top_down_unknown_values_;
std::vector<double> top_down_outlet_pressure_;
std::vector<double> top_down_child_pressure_;
```

These buffers are sized during symbolic setup from the largest top-down group.

### Replaced Scalar Top-Down Element Loop

Top-down recovery now runs through `recover_top_down_group(...)` for each top-down layer group.

The grouped recovery performs these phases:

```text
gather inlet-pressure corrections for all group lanes
write inlet-pressure delta entries
for each unknown slot:
  compute correction values for all group lanes
  scatter correction values to delta
for non-leaf groups:
  compute outlet pressure for all group lanes
  for each child slot:
    compute child inlet pressure for all group lanes
    stamp/write child inlet pressures
```

The correction formula is unchanged:

```cpp
value = slope * inlet_pressure + intercept;
```

### Preserved Stamp-Based Safety

The inlet-pressure stamp checks from Step 4 are still used. A group lane cannot read an inlet-pressure correction unless it was written during the current solve stamp.

Child inlet-pressure writes still use the same duplicate-write assertion.

## Behavior Preserved

- The top-down mathematical recovery formula is unchanged.
- Top-down layer dependency order is unchanged.
- Elements are only processed in the existing same-layer groups.
- Delta scatter still uses `replace_global_value(...)`.
- Bottom-up batch behavior from Step 7 is unchanged.
- Coefficient access behavior from Step 6 is unchanged.
- No SIMD intrinsics were added.
- The distributed tree solver was not changed.

## Guide Status

The serial vectorization guide phases are now implemented in the serial solver path:

- Phase 1: baseline measurements documented in `Optimization_Step_1.md`.
- Phase 2: serial plan/workspace SoA layout.
- Phase 3: bottom-up and top-down layer groups.
- Phase 4: avoidable serial hot-path clears removed.
- Phase 5: `2x2` and `3x3` batch dense solvers with scalar fallback.
- Phase 6: virtual coefficient-provider dispatch removed from serial hot loops.
- Phase 7: bottom-up traversal uses group assemble, batch solve, and group relation write phases.
- Phase 8: top-down traversal uses grouped lane-wise recovery phases.

The remaining structured-coefficient row-entry search in `TreeLinearization::value(...)` is still a follow-up optimization. The guide listed precomputed coefficient locations and direct tree block storage as better follow-ups rather than the minimum serial vectorization milestone.

## Intended Benefit

This step completes the grouped serial tree traversal design: both bottom-up condensation and top-down recovery now operate by layer groups instead of fully scalar per-element traversal. The loops are still explicit C++ loops, but their structure is regular and same-shape, making them more suitable for compiler vectorization and later targeted micro-optimizations.

## Verification

Verification passed:

```text
git diff --check
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

No benchmark rerun was performed in this step.
