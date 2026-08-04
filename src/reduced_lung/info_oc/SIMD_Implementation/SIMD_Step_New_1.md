# SIMD Step New 1

Date: 2026-08-04

## Goal

Add explicit batch/SIMD path accounting for the serial `TreeNewtonLinearSolver` without changing solver behavior or numerical results.

This step only adds profile counters and increments them at existing dispatch points. It does not add new SIMD kernels, change matrix assembly, change dense solve formulas, or remove scalar fallbacks.

## Files Changed

```text
src/reduced_lung/src/4C_reduced_lung_solver_profile.hpp
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
src/reduced_lung/benchmark_tests/4C_reduced_lung_solver_benchmark.cpp
```

## Profile Counters Added

Added these fields to `TreeNewtonLinearSolverProfile`:

```cpp
std::uint64_t simd_group_count = 0;
std::uint64_t simd_lane_count = 0;
std::uint64_t scalar_group_count = 0;
std::uint64_t scalar_tail_lane_count = 0;
std::uint64_t dense_fallback_count = 0;
std::uint64_t unsupported_block_fallback_count = 0;
```

Counter meanings:

- `simd_group_count`: number of batched dense-solve groups with at least one full SIMD chunk.
- `simd_lane_count`: number of element lanes processed by full SIMD chunks in batched dense solves.
- `scalar_group_count`: number of scalar dispatch points used by the serial tree solver.
- `scalar_tail_lane_count`: number of lanes processed by scalar tail loops after full SIMD chunks.
- `dense_fallback_count`: number of lanes that fell back to scalar `solve_dense_system(...)` from the batched dense-solve kernels.
- `unsupported_block_fallback_count`: number of elements solved by the generic scalar dense-solve path because their block size was not handled by the `2x2` or `3x3` batch kernels.

## Solver Instrumentation

Updated the serial tree solver batch kernels in `4C_reduced_lung_tree_linear_solver.cpp`:

```text
solve_2x2_batch(...)
solve_3x3_batch(...)
```

Both functions now receive `TreeNewtonLinearSolverProfile* profile` and record:

- One SIMD group when `simd_end > 0`.
- `simd_end` SIMD lanes for full SIMD chunks.
- `group_size - lane` scalar tail lanes before scalar tail processing.
- `fallback_count` dense fallback lanes before calling scalar `solve_dense_system(...)`.

Existing call sites now pass `profile_` into both batch-solve functions.

## Scalar Dispatch Accounting

Added `scalar_group_count` increments at existing scalar dispatch points:

- Bottom-up full scalar traversal when `use_scalar_tree_solve_` is true.
- Bottom-up generic `assemble_group(...)` path.
- Bottom-up unsupported block-size scalar solve path.
- Top-down full scalar traversal when `use_scalar_tree_solve_` is true.
- Top-down small-group scalar recovery when `group_size <= top_down_scalar_group_threshold`.
- Top-down generic `recover_top_down_group(...)` path.

Added `unsupported_block_fallback_count` for the bottom-up scalar dense-solve branch used when a group block size is neither `2` nor `3`.

## Benchmark Counters Added

Updated the serial benchmark helper `set_tree_profile_counters(...)` in:

```text
src/reduced_lung/benchmark_tests/4C_reduced_lung_solver_benchmark.cpp
```

New benchmark counters:

```text
tree_simd_groups
tree_simd_lanes
tree_scalar_groups
tree_scalar_tail_lanes
tree_dense_fallbacks
tree_unsupported_fallbacks
```

The distributed benchmark was intentionally not changed in this serial-only step.

## Behavior Preserved

- No solver formulas were changed.
- No matrix or RHS assembly logic was changed.
- No tree traversal order was changed.
- No scalar fallback was removed.
- No distributed tree solver code was changed.

The new counters are only updated when a `TreeNewtonLinearSolverProfile` is attached.

## Expected Use

For a large serial structured `NewtonTree` benchmark, these counters should now make it visible whether the run is using the existing batch dense-solve kernels:

```text
tree_simd_groups > 0
tree_simd_lanes > 0
```

They also expose remaining scalar work:

```text
tree_scalar_groups
```

For current production block sizes, `tree_unsupported_fallbacks` should remain zero in normal serial `NewtonTree` runs.

## Verification

Suggested checks:

```text
git diff --check
cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
```

Also run a serial `NewtonTree` benchmark or generation-10 input and inspect the new benchmark counters.
