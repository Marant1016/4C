# SIMD Step New 8

Date: 2026-08-04

## Goal

Gate the serial `TreeNewtonLinearSolver` small-group scalar bypasses so production defaults remain unchanged while benchmarks or targeted tests can force grouped batch traversal.

Before this step, two scalar bypasses were unconditional:

- Whole-tree scalar traversal for trees with at most seven elements.
- Scalar top-down recovery for groups with at most two lanes.

## Files Changed

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
```

## New Context Flag

Added a disabled-by-default serial solver context option:

```cpp
bool force_batch_tree_solve = false;
```

This is stored on `TreeNewtonLinearSolver` as:

```cpp
bool force_batch_tree_solve_ = false;
```

The default keeps existing small-tree and tiny top-down scalar behavior unchanged.

## Whole-Tree Scalar Gate

The small-tree scalar solve selection now honors the force-batch flag:

```cpp
use_scalar_tree_solve_ =
    !force_batch_tree_solve_ && element_count <= scalar_tree_element_threshold;
```

When `force_batch_tree_solve` is `true`, small trees use the grouped bottom-up and top-down paths instead of the element-by-element scalar traversal.

## Top-Down Small-Group Gate

The tiny top-down group bypass now also honors the force-batch flag:

```cpp
const bool use_scalar_top_down_group =
    !force_batch_tree_solve_ && group_size <= top_down_scalar_group_threshold;
```

When `force_batch_tree_solve` is `true`, supported `2x2` and `3x3` groups route through the dedicated top-down recovery paths even if the group has only one or two lanes.

## Behavior Preserved

- Existing solver construction keeps `force_batch_tree_solve = false` by default.
- Production small-tree scalar traversal remains unchanged unless the new flag is explicitly set.
- Unsupported future block sizes still use scalar fallback paths.
- Sparse-Jacobian validation remains available.
- The distributed tree solver was not changed.

## Intended Use

The new flag is intended for benchmark or targeted validation runs that need strict grouped traversal and meaningful batch-path accounting. In normal production runs, the existing scalar shortcuts remain active for tiny cases.

## Verification

Commands run after implementation:

```text
cmake --build build/debug --target unittests_reduced_lung --parallel 4
git diff --check
ctest -R "^unittests_reduced_lung$" --output-on-failure
```
