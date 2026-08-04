# SIMD Step New 3

Date: 2026-08-04

## Goal

Convert the serial `TreeNewtonLinearSolver` `2x2` batch-solve path to use direct structure-of-arrays (SoA) workspaces instead of reading the interleaved per-element `workspace_matrix_`, `workspace_rhs_constant_`, and `workspace_rhs_inlet_pressure_` arrays inside `solve_2x2_batch(...)`.

This step keeps downstream solver behavior unchanged by copying the `2x2` SoA solve results back into the existing `workspace_intercept_` and `workspace_slope_` arrays after the batch solve.

## Files Changed

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
```

## New 2x2 SoA Workspaces

Added these members to `TreeNewtonLinearSolver`:

```cpp
std::vector<double> batch_2x2_a00_;
std::vector<double> batch_2x2_a01_;
std::vector<double> batch_2x2_a10_;
std::vector<double> batch_2x2_a11_;
std::vector<double> batch_2x2_rhs_constant0_;
std::vector<double> batch_2x2_rhs_constant1_;
std::vector<double> batch_2x2_rhs_inlet_pressure0_;
std::vector<double> batch_2x2_rhs_inlet_pressure1_;
std::vector<double> batch_2x2_intercept0_;
std::vector<double> batch_2x2_intercept1_;
std::vector<double> batch_2x2_slope0_;
std::vector<double> batch_2x2_slope1_;
```

They are allocated in `build_symbolic_plan()` using `max_2x2_group_size`, alongside the existing `batch_2x2_fallback_lanes_` workspace.

## Batch Solve Changes

Updated `solve_2x2_batch(...)` so it now consumes and writes SoA arrays directly:

```text
batch_2x2_a00_ ... batch_2x2_a11_
batch_2x2_rhs_constant0_ ... batch_2x2_rhs_constant1_
batch_2x2_rhs_inlet_pressure0_ ... batch_2x2_rhs_inlet_pressure1_
batch_2x2_intercept0_ ... batch_2x2_intercept1_
batch_2x2_slope0_ ... batch_2x2_slope1_
```

The SIMD loop now gathers by lane index from these SoA arrays rather than by grouped element offsets into interleaved workspaces.

The scalar fallback inside `solve_2x2_batch(...)` now builds a small stack-local `2x2` dense system from the SoA lane values, calls `solve_dense_system(...)`, and writes the fallback result back into the SoA output arrays.

## Structured 2x2 Assembly Changes

The structured SIMD branches for these existing grouped assembly paths now write directly into the `2x2` SoA workspaces:

```text
assemble_2x2_leaf_group(...)
assemble_2x2_one_child_group(...)
assemble_2x2_two_child_group(...)
```

Changes inside the structured chunk helpers:

- Equation RHS values write to `batch_2x2_rhs_constant0_` and `batch_2x2_rhs_constant1_`.
- Inlet-pressure RHS coefficients write to `batch_2x2_rhs_inlet_pressure0_` and `batch_2x2_rhs_inlet_pressure1_`.
- Matrix entries write to `batch_2x2_a00_`, `batch_2x2_a01_`, `batch_2x2_a10_`, and `batch_2x2_a11_`.
- Child downstream-flow matrix updates now add directly to `batch_2x2_a10_` or `batch_2x2_a11_`, depending on the parent outlet-pressure unknown index.
- Child downstream-flow RHS shifts subtract directly from `batch_2x2_rhs_constant1_`.

## Compatibility For Existing Scalar Paths

Added `pack_2x2_group_from_workspace(...)` for paths that still assemble `2x2` blocks through the interleaved workspaces.

This pack path is used when direct structured SoA assembly is not active, including:

- Sparse-Jacobian coefficient source.
- Builds without `<experimental/simd>`.
- Any `2x2` path that still relies on scalar workspace assembly.

This preserves existing validation behavior while allowing the same SoA solve kernel to be used afterward.

## Result Write-Back

Added `write_2x2_batch_solution_to_workspace(...)` after `solve_2x2_batch(...)`.

It copies:

```text
batch_2x2_intercept0_ -> workspace_intercept_[unknown_begin]
batch_2x2_intercept1_ -> workspace_intercept_[unknown_begin + 1]
batch_2x2_slope0_     -> workspace_slope_[unknown_begin]
batch_2x2_slope1_     -> workspace_slope_[unknown_begin + 1]
```

This keeps downstream code unchanged for now:

- `write_subtree_relation_group(...)`
- top-down recovery paths
- generic fallback paths

Future steps can remove this copy by making subtree relation writes and top-down recovery read the `2x2` SoA outputs directly.

## Behavior Preserved

- No mathematical solve formulas changed.
- Existing dense fallback behavior remains available for real lanes.
- Existing scalar workspace assembly remains available for sparse and non-SIMD paths.
- Existing subtree relation and top-down recovery code still read `workspace_intercept_` and `workspace_slope_`.
- The distributed tree solver was not changed.

## Expected Performance Impact

For structured serial `2x2` groups, the hot batch solve no longer gathers matrix/RHS values indirectly through per-element `matrix_begin` and `unknown_begin` offsets.

Expected benchmark signal:

```text
tree_dense_s should improve for 2x2-heavy serial NewtonTree cases.
tree_simd_lanes should remain comparable to Step 2.
tree_scalar_tail_lanes should remain zero for SIMD-enabled 2x2 batch solves.
```

## Verification

Commands run after this step:

```text
cmake --build build/debug --target unittests_reduced_lung --parallel 4
```

Recommended final checks:

```text
git diff --check
ctest -R "^unittests_reduced_lung$" --output-on-failure
```
