# SIMD Step New 4

Date: 2026-08-04

## Goal

Add structured SIMD batch assembly for serial `TreeNewtonLinearSolver` `3x3` groups so structured large-tree runs no longer assemble current production `3x3` element groups through the generic scalar `assemble_group(...)` path.

This step covers:

- `3x3 + 0 child`
- `3x3 + 1 child`
- `3x3 + 2 children`

The sparse-Jacobian coefficient source and non-SIMD builds keep the existing scalar workspace assembly behavior.

## Files Changed

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
```

## Structured 3x3 Batch Assembly

Added structured SIMD chunk helpers inside `TreeNewtonLinearSolver::solve(...)`:

```text
assemble_3x3_equation_rows_structured_chunk(...)
add_3x3_child_contribution_structured_chunk(...)
subtract_3x3_rhs_shift_structured_chunk(...)
```

The equation-row helper gathers three RHS values, three inlet-pressure coefficients, and nine local matrix coefficients per SIMD chunk. It writes directly into the existing `batch_3x3_*` SoA arrays.

The child-contribution helper mirrors the existing `2x2` structured path. It computes the child pressure relation, stores `child_pressure_slope_` and `child_pressure_intercept_`, uses the child subtree relation, and adds the downstream flow contribution to the last `3x3` matrix row.

## Group Dispatch

Added grouped assembly wrappers:

```text
assemble_3x3_leaf_group(...)
assemble_3x3_one_child_group(...)
assemble_3x3_two_child_group(...)
```

The bottom-up grouped dispatch now routes structured `3x3` production shapes to these wrappers before the generic scalar fallback.

## Direct SoA Solve Guard

Updated `solve_3x3_batch(...)` with an `inputs_already_packed` flag.

When structured SIMD assembly fills `batch_3x3_*` directly, the solve skips the old scalar copy from `workspace_matrix_`, `workspace_rhs_constant_`, and `workspace_rhs_inlet_pressure_` into SoA arrays. This is required because the direct structured path does not populate those interleaved workspace inputs.

The non-structured path still uses the old scalar pack from workspace before the SIMD `3x3` solve.

## Dense Fallback

The `3x3` dense fallback now builds a stack-local dense system from the SoA lane values, solves it with `solve_dense_system(...)`, and writes results through the existing `write_solution(...)` path.

This keeps fallback support available for both direct structured SoA assembly and sparse/workspace assembly.

## Behavior Preserved

- Numerical solve formulas for the direct `3x3` batch path are unchanged.
- Sparse-Jacobian validation remains available through scalar workspace assembly.
- Builds without `<experimental/simd>` still use scalar workspace assembly for `3x3` groups.
- `workspace_intercept_` and `workspace_slope_` are still populated after `3x3` solves, preserving downstream subtree relation and top-down recovery code.
- The distributed tree solver was not changed.

## Verification

Commands run after implementation:

```text
cmake --build build/debug --target unittests_reduced_lung --parallel 4
git diff --check
ctest -R "^unittests_reduced_lung$" --output-on-failure
```
