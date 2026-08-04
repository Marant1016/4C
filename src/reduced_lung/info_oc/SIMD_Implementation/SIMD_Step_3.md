# SIMD Step 3

Date: 2026-08-03

## Goal

Add an explicit SIMD fast path to the serial `TreeNewtonLinearSolver` `2x2` dense batch solve.

This is the first step that uses the SIMD helper layer from Step 1 and the grouped offset caches from Step 2 in an actual numeric kernel.

## Files Changed

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
```

## Implementation

Updated `solve_2x2_batch(...)` to accept grouped offset caches:

```cpp
const std::vector<int>& grouped_unknown_begin
const std::vector<int>& grouped_matrix_begin
```

The function no longer recomputes `unknown_begin` and `matrix_begin` through:

```cpp
element_index -> unknown_offset_[element_index]
element_index -> matrix_offset_[element_index]
```

Instead, it reads directly from the grouped caches:

```cpp
const int unknown_begin = grouped_unknown_begin[grouped_index];
const int matrix_begin = grouped_matrix_begin[grouped_index];
```

The bottom-up `2x2` call site now passes:

```cpp
grouped_unknown_begin_
grouped_matrix_begin_
```

## SIMD Fast Path

When `<experimental/simd>` is available, `solve_2x2_batch(...)` processes full chunks of:

```cpp
tree_solver_simd::width()
```

For each SIMD chunk, it gathers:

```text
a00, a01, a10, a11
rhs_constant0, rhs_constant1
rhs_inlet_pressure0, rhs_inlet_pressure1
```

Then it computes the direct `2x2` solution formulas with SIMD arithmetic:

```text
det = a00 * a11 - a01 * a10
inv_det = 1.0 / det

intercept0 = (rhs_constant0 * a11 - a01 * rhs_constant1) * inv_det
intercept1 = (a00 * rhs_constant1 - rhs_constant0 * a10) * inv_det

slope0 = (rhs_inlet_pressure0 * a11 - a01 * rhs_inlet_pressure1) * inv_det
slope1 = (a00 * rhs_inlet_pressure1 - rhs_inlet_pressure0 * a10) * inv_det
```

The SIMD results are scattered back to:

```text
workspace_intercept_
workspace_slope_
```

through the existing `intercept` and `slope` vectors passed to the batch solver.

## Fallback Behavior

The SIMD path uses a conservative determinant check:

```text
all lanes must be finite
all lanes must satisfy abs(det) > pivot_tolerance
```

If any lane in a SIMD chunk fails this check, the entire chunk is sent to the existing pivoted scalar fallback:

```cpp
solve_dense_system(..., block_size = 2, ...)
```

This is intentionally conservative. A later optimization can fallback only invalid lanes, but full-chunk fallback keeps this first SIMD numeric step simple and safe.

Scalar fallback remains active for:

- unsafe SIMD chunks,
- tail elements after full SIMD chunks,
- group sizes smaller than SIMD width,
- builds where `<experimental/simd>` is unavailable.

## Behavior Preserved

- The small-tree scalar path is unchanged.
- Group construction and traversal order are unchanged.
- `2x2` mathematical formulas are unchanged for valid pivots.
- Unsafe `2x2` systems still use the existing pivoted dense solver.
- `3x3` batch solve is unchanged.
- Bottom-up assembly is unchanged.
- Top-down recovery is unchanged.
- The distributed tree solver is unchanged.

## Expected Benefit

This step targets `tree_dense_s` for large grouped `2x2` cases.

The largest benefit is expected when a bottom-up layer group contains many `2x2` elements and most determinants pass the fast-path check.

Small models should remain scalar through the existing small-tree path or scalar tails, so this step should not meaningfully affect tiny cases.

## Verification

Verification target for this step:

```text
git diff --check
cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

Performance check target after a release build is available:

```text
./build/release/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/LinearSolve/BalancedAirways/StructuredTree
```

Record at least:

```text
tree_solve_s
tree_bottom_up_s
tree_dense_s
tree_dense_solves
```

## Follow-Up

The next likely step is SIMD-aware structured `2x2` assembly for:

```text
assemble_2x2_leaf_group(...)
assemble_2x2_one_child_group(...)
assemble_2x2_two_child_group(...)
```

That should improve the non-dense-solve part of `tree_bottom_up_s`.
