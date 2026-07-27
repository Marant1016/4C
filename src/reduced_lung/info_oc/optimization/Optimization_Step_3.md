# Optimization Step 3 - Serial Tree Solver Layer Groups

## Scope

This step implements Phase 3 from the serial tree solver vectorization guide.

Only the serial `TreeNewtonLinearSolver` traversal data was changed. The distributed `DistributedTreeNewtonLinearSolver` implementation was intentionally left unchanged.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Added Element Groups

The serial solver now stores same-shape element groups for bottom-up and top-down traversal:

```cpp
struct ElementGroup
{
  int begin = 0;
  int end = 0;
  int block_size = 0;
  int child_count = 0;
};

std::vector<int> grouped_element_indices_;
std::vector<std::vector<ElementGroup>> bottom_up_layer_groups_;
std::vector<std::vector<ElementGroup>> top_down_layer_groups_;
```

Each group range indexes into `grouped_element_indices_`.

### Built Groups During Symbolic Setup

After the Step 2 SoA arrays are populated, symbolic setup now builds grouped views of the existing tree layers.

Elements are grouped within each layer by:

```text
block_size_[element_index]
child_interface_count_[element_index]
```

Shape keys are sorted deterministically by `block_size`, then `child_count`.

The existing tree layers are still the source of truth. This step only creates a grouped traversal view inside each layer.

### Switched Serial Traversal To Grouped Ranges

The bottom-up and top-down solve phases now iterate through grouped layer data:

```cpp
for (const auto& layer_groups : bottom_up_layer_groups_)
{
  for (const auto& group : layer_groups)
  {
    for (int grouped_index = group.begin; grouped_index < group.end; ++grouped_index)
    {
      const int element_index = grouped_element_indices_[grouped_index];
      // Existing scalar per-element work.
    }
  }
}
```

The per-element solve body remains scalar and mathematically unchanged.

## Behavior Preserved

- The original bottom-up and top-down layer ordering is preserved.
- Elements are only reordered within a single layer.
- The scalar local solve and substitution algorithm is unchanged.
- The equation order and unknown order are unchanged.
- The coefficient-provider paths are unchanged.
- No vectorized dense kernels were added.
- No direct structured coefficient lookup path was added.
- Existing array clears were intentionally kept for later phases.
- The distributed tree solver was not changed.

## Intended Benefit

This step prepares the serial solver for later batch processing. Same-shape elements are now contiguous in traversal order, which allows later phases to gather coefficients, substitute child relations, and call shape-specific batch kernels without first building grouping data during the hot solve path.

This step is not expected to provide the main speedup by itself. It isolates traversal grouping from later numerical-kernel changes.

## Verification

Verification passed:

```text
git diff --check
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

No benchmark rerun was performed in this step.
