# SIMD Step 2

Date: 2026-08-03

## Goal

Add grouped offset caches parallel to `grouped_element_indices_` for the serial `TreeNewtonLinearSolver`.

This prepares the grouped traversal for later SIMD kernels by avoiding repeated per-lane lookups from `element_index` to hot offset arrays.

## Files Changed

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
```

## Implementation

Added grouped offset cache members to `TreeNewtonLinearSolver`:

```cpp
std::vector<int> grouped_unknown_begin_;
std::vector<int> grouped_matrix_begin_;
std::vector<int> grouped_equation_begin_;
std::vector<int> grouped_child_begin_;
```

These arrays are parallel to:

```cpp
std::vector<int> grouped_element_indices_;
```

For every grouped element entry, the caches store:

```text
grouped_unknown_begin_[i]   == unknown_offset_[element_index]
grouped_matrix_begin_[i]    == matrix_offset_[element_index]
grouped_equation_begin_[i]  == equation_offset_[element_index]
grouped_child_begin_[i]     == child_interface_offset_[element_index]
```

The caches are cleared with `grouped_element_indices_` in `build_symbolic_plan()`.

They are populated inside `build_layer_groups(...)` immediately after appending an element to `grouped_element_indices_`:

```cpp
grouped_element_indices_.push_back(element_index);
grouped_unknown_begin_.push_back(unknown_offset_[element_index_size]);
grouped_matrix_begin_.push_back(matrix_offset_[element_index_size]);
grouped_equation_begin_.push_back(equation_offset_[element_index_size]);
grouped_child_begin_.push_back(child_interface_offset_[element_index_size]);
```

## Validation

Added validation that all grouped cache arrays have the same size as `grouped_element_indices_`.

Extended `validate_grouped_traversal(...)` to check each cached offset against the corresponding element plan entry.

This validates both grouped traversals:

```text
bottom_up_layer_groups_
top_down_layer_groups_
```

## Behavior Preserved

- The small-tree scalar path is unchanged.
- The large-tree grouped traversal order is unchanged.
- No solver kernel reads the new caches yet.
- `solve_2x2_batch(...)` is unchanged.
- `solve_3x3_batch(...)` is unchanged.
- `assemble_2x2_*_group(...)` functions are unchanged.
- `recover_2x2_top_down_group(...)` is unchanged.
- The distributed tree solver is unchanged.
- No explicit SIMD arithmetic was added in this step.

This step should produce no numerical or performance behavior change by itself.

## Purpose For Next Steps

Later SIMD kernels can use these caches to avoid repeated setup like:

```cpp
const int element_index = grouped_element_indices_[grouped_index];
const int unknown_begin = unknown_offset_[element_index];
const int matrix_begin = matrix_offset_[element_index];
const int equation_begin = equation_offset_[element_index];
const int child_begin = child_interface_offset_[element_index];
```

Instead, they can load directly from the grouped cache arrays:

```cpp
const int unknown_begin = grouped_unknown_begin_[grouped_index];
const int matrix_begin = grouped_matrix_begin_[grouped_index];
const int equation_begin = grouped_equation_begin_[grouped_index];
const int child_begin = grouped_child_begin_[grouped_index];
```

The first intended consumer remains the structured `2x2` grouped path, starting with SIMD support in `solve_2x2_batch(...)`.

## Verification

Verification target for this step:

```text
git diff --check
cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

Because this step only adds unused grouped offset caches and validation, correctness should be identical to the previous implementation.
