# SIMD Step New 9

Date: 2026-08-04

## Goal

Define the numerical fallback policy for serial `TreeNewtonLinearSolver` batched dense solves.

The policy keeps scalar dense fallback as a correctness path by default, counts fallback lanes through profiling, and adds an opt-in strict mode for benchmark or validation runs that must prove they do not rely on scalar fallback.

## Files Changed

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
```

## New Context Flag

Added a disabled-by-default serial solver context option:

```cpp
bool error_on_dense_fallback = false;
```

This is stored on `TreeNewtonLinearSolver` as:

```cpp
bool error_on_dense_fallback_ = false;
```

Default behavior is unchanged: dense fallback remains allowed and counted.

## Batch Solver Policy

Both `solve_2x2_batch(...)` and `solve_3x3_batch(...)` now receive the strict fallback flag.

The existing fallback detection remains unchanged:

- `2x2` fallback is requested for non-finite or near-singular determinants.
- `3x3` fallback is requested for unsafe pivots, non-finite factorization values, or non-finite solutions.

The existing profile accounting remains in place:

```cpp
profile->dense_fallback_count += fallback_count;
```

If `error_on_dense_fallback` is enabled and `fallback_count > 0`, the solver now throws before running scalar dense fallback.

## Default Fallback Behavior

When `error_on_dense_fallback` is false, fallback lanes still use stack-local dense systems and `solve_dense_system(...)`, preserving the existing correctness behavior.

This applies to both `2x2` and `3x3` batch solve paths.

## Strict Mode Use

Strict mode is intended for benchmark or targeted validation runs where scalar fallback should be treated as a numerical or performance failure.

Recommended interpretation:

- `dense_fallback_count == 0` is expected for normal well-conditioned validation and benchmark cases.
- `dense_fallback_count > 0` in normal big-tree cases should be investigated.
- Batched pivoting remains deferred; scalar fallback is still the default correctness policy.

## Behavior Preserved

- Production default behavior is unchanged.
- Dense fallback remains available unless strict mode is explicitly enabled.
- Existing fallback counters remain visible in profiles and benchmark counters.
- Sparse-Jacobian validation remains available.
- The distributed tree solver was not changed.

## Verification

Commands run after implementation:

```text
cmake --build build/debug --target unittests_reduced_lung --parallel 4
git diff --check
ctest -R "^unittests_reduced_lung$" --output-on-failure
```
