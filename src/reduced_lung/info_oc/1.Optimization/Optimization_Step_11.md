# Optimization Step 11 - Reduce Top-Down Staging Overhead

## Scope

This step implements the top-down staging overhead follow-up from the 2026-07-31 serial tree solver guide.

Only the serial `TreeNewtonLinearSolver` top-down recovery path was changed. The distributed `DistributedTreeNewtonLinearSolver` implementation was intentionally left unchanged.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Added Scalar Top-Down Element Recovery

The serial solver now has a scalar `recover_top_down_element(...)` path for small top-down groups.

The scalar path performs the same recovery operations directly for one element:

```text
read inlet-pressure correction
write inlet-pressure correction to delta
compute and write each local unknown correction
if the element has children:
  compute outlet-pressure correction
  compute and stamp each child inlet-pressure correction
```

It avoids the grouped path's temporary staging buffers:

```cpp
top_down_inlet_pressure_
top_down_unknown_values_
top_down_outlet_pressure_
top_down_child_pressure_
```

### Added Small-Group Dispatch Threshold

Top-down recovery now dispatches groups with size `<= 2` to scalar element recovery:

```cpp
constexpr int top_down_scalar_group_threshold = 2;
```

Larger groups continue to use the existing grouped `recover_top_down_group(...)` path.

## Behavior Preserved

- The top-down mathematical recovery formula is unchanged.
- Top-down layer dependency order is unchanged.
- The grouped path remains active for larger groups.
- Inlet-pressure stamp validation is unchanged.
- Duplicate child inlet-pressure write assertions are unchanged.
- `delta.replace_global_value(...)` remains the correction scatter operation.
- Bottom-up condensation and coefficient access behavior from Steps 9-10 are unchanged.
- The distributed tree solver was not changed.

## Intended Benefit

Step 8 made top-down recovery vectorization-friendly by staging group lanes through reusable buffers. The Step 8 benchmark rerun showed that this staging can be too expensive for tiny groups.

This step keeps the grouped path where it can still amortize staging overhead, but restores a low-overhead direct path for singleton and two-lane groups. The expected benefit is lower `tree_top_down_s` for small and medium structured trees without changing larger-group behavior.

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

The next release benchmark run should compare `tree_top_down_s` against Step 10 and Step 8, especially for `StructuredTree/2_mean`, `StructuredTree/3_mean`, and `SingleTerminalUnit/NewtonTree_mean`.
