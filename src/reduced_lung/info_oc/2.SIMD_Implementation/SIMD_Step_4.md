# SIMD Step 4

Date: 2026-08-03

## Goal

Add SIMD-aware structured assembly for serial `TreeNewtonLinearSolver` `2x2` bottom-up groups.

This step extends SIMD use beyond the `2x2` dense solve from Step 3 and targets the assembly phase immediately before `solve_2x2_batch(...)`.

## Files Changed

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
```

## Implementation

Added guarded SIMD helpers inside the serial tree solve path when `<experimental/simd>` is available:

```text
require_structured_coefficients(...)
assemble_2x2_equation_rows_structured_chunk(...)
add_2x2_child_contribution_structured_chunk(...)
subtract_2x2_rhs_shift_structured_chunk(...)
```

These helpers are only used when:

```cpp
coefficient_source_ == TreeNewtonLinearSolverCoefficientSource::StructuredTreeBlocks
```

The sparse-Jacobian path stays on the existing scalar assembly path.

## SIMD Equation-Row Assembly

`assemble_2x2_equation_rows_structured_chunk(...)` processes one full SIMD chunk of grouped `2x2` elements.

It gathers:

```text
rhs row 0
rhs row 1
inlet-pressure coefficient row 0
inlet-pressure coefficient row 1
matrix coefficients a00, a01, a10, a11
```

Then it scatters assembled values into:

```text
workspace_rhs_constant_
workspace_rhs_inlet_pressure_
workspace_matrix_
```

It uses grouped caches from Step 2:

```text
grouped_unknown_begin_
grouped_matrix_begin_
grouped_equation_begin_
```

This avoids recomputing element offsets through `grouped_element_indices_` for each lane.

## SIMD Child Contribution Assembly

`add_2x2_child_contribution_structured_chunk(...)` vectorizes the child-substitution arithmetic for `2x2` groups with one or two children.

For each child slot, it computes across a SIMD chunk:

```text
child_pressure_slope = -pressure_parent_coeff / pressure_child_coeff
child_pressure_intercept = pressure_rhs / pressure_child_coeff

child_flow_slope = subtree_relation_G * child_pressure_slope
child_flow_intercept = subtree_relation_G * child_pressure_intercept + subtree_relation_h

matrix_update = flow_child_coeff * child_flow_slope
rhs_update = flow_child_coeff * child_flow_intercept
```

Then it scatters:

```text
child_pressure_slope_
child_pressure_intercept_
workspace_matrix_
```

and accumulates a vector `rhs_shift` for the chunk.

`subtract_2x2_rhs_shift_structured_chunk(...)` subtracts that vector shift from the second `2x2` RHS entry for each grouped element.

## Dispatch

The existing group dispatch shape is preserved:

```text
assemble_2x2_leaf_group(...)
assemble_2x2_one_child_group(...)
assemble_2x2_two_child_group(...)
```

Each helper now has this structure:

```text
if structured coefficient source and SIMD is available:
  SIMD full chunks
  scalar tail
  return

existing scalar assembly path
```

## Fallback Behavior

Scalar assembly remains active for:

- small trees using `use_scalar_tree_solve_`,
- sparse-Jacobian coefficient source,
- builds without `<experimental/simd>`,
- tail elements after full SIMD chunks,
- non-`2x2` groups.

Required structured child coefficients are still validated. If any lane in a chunk has a missing or near-zero required coefficient, the solver raises the same kind of error as the scalar path.

## Behavior Preserved

- The small-tree scalar path is unchanged.
- The sparse-Jacobian path is unchanged.
- The grouped traversal order is unchanged.
- `solve_2x2_batch(...)` behavior from Step 3 is unchanged.
- `3x3` assembly and solve paths are unchanged.
- Top-down recovery is unchanged.
- The distributed tree solver is unchanged.

## Expected Benefit

This step targets the assembly portion of `tree_bottom_up_s` for large structured `2x2` groups.

The best expected improvement is for large groups with `child_count == 1` or `child_count == 2`, because those paths include child pressure/flow substitution arithmetic in addition to equation-row loading.

Leaf groups benefit mainly from batched coefficient/RHS loading and stores.

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
tree_top_down_s
```

## Follow-Up

The next likely SIMD step is `2x2` top-down recovery:

```text
recover_2x2_top_down_group(...)
```

That should target `tree_top_down_s` by vectorizing correction recovery and child inlet-pressure propagation for structured `2x2` groups.
