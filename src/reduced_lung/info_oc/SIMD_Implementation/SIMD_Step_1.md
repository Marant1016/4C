# SIMD Step 1

Date: 2026-08-03

## Goal

Add an isolated SIMD helper layer for the serial `TreeNewtonLinearSolver` without changing solver behavior.

This is preparation work for later SIMD kernels. No bottom-up assembly, dense solve, top-down recovery, dispatch threshold, or scalar fallback behavior was changed in this step.

## Files Changed

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
```

## Implementation

Added a guarded include for C++ experimental SIMD:

```cpp
#if defined(__has_include)
#if __has_include(<experimental/simd>)
#include <experimental/simd>
#define FOUR_C_REDUCED_LUNG_HAS_EXPERIMENTAL_SIMD 1
#else
#define FOUR_C_REDUCED_LUNG_HAS_EXPERIMENTAL_SIMD 0
#endif
#else
#define FOUR_C_REDUCED_LUNG_HAS_EXPERIMENTAL_SIMD 0
#endif
```

Added an internal helper namespace inside the anonymous namespace of `4C_reduced_lung_tree_linear_solver.cpp`:

```cpp
namespace tree_solver_simd
{
  ...
}
```

When `<experimental/simd>` is available, the helper provides:

```text
tree_solver_simd::Double
tree_solver_simd::Mask
tree_solver_simd::available
tree_solver_simd::width()
tree_solver_simd::full_chunk_end(...)
tree_solver_simd::has_full_chunk(...)
tree_solver_simd::gather(...)
tree_solver_simd::scatter(...)
```

When `<experimental/simd>` is not available, the helper provides a scalar-compatible fallback state:

```text
tree_solver_simd::available == false
tree_solver_simd::width() == 1
tree_solver_simd::full_chunk_end(...) returns begin
tree_solver_simd::has_full_chunk(...) returns false
```

The fallback path intentionally does not expose SIMD vector types. Future SIMD kernels should guard SIMD-only code with the compile-time feature macro or `tree_solver_simd::available`.

## Behavior Preserved

- The current small-tree scalar threshold is unchanged.
- The current grouped large-tree traversal is unchanged.
- `solve_2x2_batch(...)` is unchanged.
- `solve_3x3_batch(...)` is unchanged.
- `assemble_2x2_*_group(...)` functions are unchanged.
- `recover_2x2_top_down_group(...)` is unchanged.
- The distributed tree solver is unchanged.
- No explicit SIMD kernel is called yet.

This step should produce no numerical or performance behavior change by itself.

## Purpose For Next Steps

Later steps can use this helper to process full SIMD chunks inside same-shape element groups:

```text
for each ElementGroup:
  SIMD chunks of tree_solver_simd::width() elements
  scalar tail for remaining elements
```

The first intended consumer is the structured `2x2` grouped path, starting with `solve_2x2_batch(...)`.

## Verification

Verification target for this step:

```text
git diff --check
cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

Because this step only adds unused helper scaffolding, correctness should be identical to the previous grouped scalar implementation.
