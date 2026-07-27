# Optimization Step 2 - Serial Tree Solver SoA Layout

## Scope

This step implements Phase 2 from the serial tree solver vectorization guide.

Only the serial `TreeNewtonLinearSolver` path was changed. The distributed `DistributedTreeNewtonLinearSolver` implementation was intentionally left unchanged.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Converted Serial Plan Data To SoA

The serial `ElementSolvePlan` storage was removed from `TreeNewtonLinearSolver` and replaced by fixed per-element arrays:

```cpp
std::vector<int> global_element_id_;
std::vector<int> inlet_pressure_local_dof_;
std::vector<int> inlet_flow_unknown_index_;
std::vector<int> outlet_pressure_unknown_index_;
std::vector<int> block_size_;
std::vector<int> child_interface_count_;
std::vector<unsigned char> is_leaf_;
```

Variable-size per-element data now uses flat arrays with offsets:

```cpp
std::vector<int> unknown_offset_;
std::vector<int> equation_offset_;
std::vector<int> child_interface_offset_;
std::vector<int> matrix_offset_;

std::vector<int> unknown_global_dof_ids_;
std::vector<int> unknown_local_dof_ids_;
std::vector<int> equation_rows_;
```

Child-interface metadata was also flattened:

```cpp
std::vector<int> child_element_index_;
std::vector<int> pressure_row_;
std::vector<int> parent_outlet_pressure_local_dof_;
std::vector<int> child_inlet_pressure_local_dof_;
std::vector<int> child_inlet_flow_local_dof_;
std::vector<int> parent_outlet_pressure_unknown_index_;
```

### Converted Serial Workspace Data To Flat Buffers

The serial `ElementWorkspace` storage was removed and replaced by contiguous workspace arrays:

```cpp
std::vector<double> workspace_matrix_;
std::vector<double> workspace_rhs_constant_;
std::vector<double> workspace_rhs_inlet_pressure_;
std::vector<double> workspace_intercept_;
std::vector<double> workspace_slope_;
std::vector<double> child_pressure_slope_;
std::vector<double> child_pressure_intercept_;
```

Each element accesses its workspace through the offsets built during symbolic setup.

### Kept Debug Context Out Of Hot Plan Data

The old per-plan `context` string was moved to a separate array:

```cpp
std::vector<std::string> element_context_;
```

This keeps diagnostic text separate from the hot numeric plan arrays.

### Updated Scalar Solve To Use SoA Offsets

The bottom-up and top-down loops still traverse `tree_metadata_.bottom_up_layers` and `tree_metadata_.top_down_layers` exactly as before, but now access plan and workspace data through offset ranges instead of per-element structs.

The dense local solver now accepts contiguous workspace slices:

```cpp
solve_dense_system(
    workspace_matrix_.data() + matrix_begin,
    workspace_rhs_constant_.data() + unknown_begin,
    workspace_rhs_inlet_pressure_.data() + unknown_begin,
    workspace_intercept_.data() + unknown_begin,
    workspace_slope_.data() + unknown_begin,
    block_size,
    pivot_tolerance_,
    element_context_[element_index]);
```

## Behavior Preserved

- The scalar bottom-up/top-down tree algorithm is unchanged.
- The equation order and unknown order are unchanged.
- The pivoted dense local solve is unchanged mathematically.
- The sparse and structured coefficient-provider paths are unchanged.
- Existing full-array clears were intentionally kept for later phases.
- No layer grouping was added in this step.
- No vectorized dense kernels or SIMD intrinsics were added.
- No direct structured coefficient lookup path was added.
- The distributed tree solver was not refactored.

## Intended Benefit

This step prepares the serial tree solver for later layer grouping and batch processing. By flattening plan and workspace data, same-shape element batches can later be represented as contiguous ranges without first extracting data from per-element structs.

This step is not expected to be the main performance win by itself. Its purpose is to isolate data-layout risk before adding grouped traversal or specialized dense kernels.

## Verification

Verification passed:

```text
git diff --check
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

No benchmark rerun was performed in this step.
