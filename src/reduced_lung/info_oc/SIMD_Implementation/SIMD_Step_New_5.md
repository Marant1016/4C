# SIMD Step New 5

Date: 2026-08-04

## Goal

Make `solve_3x3_batch(...)` a fully batch-oriented SoA kernel. After Step 4, structured `3x3` assembly could already fill the `batch_3x3_*` arrays directly, but the solve function still carried the older mixed workspace/SoA contract.

This step removes that mixed contract and makes workspace packing an explicit compatibility step outside the solve kernel.

## Files Changed

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
```

## 3x3 Batch Solve Contract

Updated `solve_3x3_batch(...)` so it now consumes only:

```text
batch_3x3_a00_ ... batch_3x3_a22_
batch_3x3_rhs_constant0_ ... batch_3x3_rhs_constant2_
batch_3x3_rhs_inlet_pressure0_ ... batch_3x3_rhs_inlet_pressure2_
batch_3x3_intercept0_ ... batch_3x3_intercept2_
batch_3x3_slope0_ ... batch_3x3_slope2_
```

Removed from the solve function contract:

- `workspace_matrix_`
- `workspace_rhs_constant_`
- `workspace_rhs_inlet_pressure_`
- `workspace_intercept_`
- `workspace_slope_`
- `grouped_unknown_begin_`
- `grouped_matrix_begin_`
- `inputs_already_packed`

The solve function no longer performs a lane-by-lane copy from interleaved workspace arrays before solving, and it no longer writes solutions directly back to `workspace_intercept_` or `workspace_slope_`.

## Workspace Compatibility Path

Added `pack_3x3_group_from_workspace(...)` for sparse-Jacobian validation and non-direct `3x3` paths.

This helper mirrors the existing `2x2` pack helper and copies interleaved workspace data into the `3x3` SoA arrays before calling `solve_3x3_batch(...)`.

Structured SIMD `3x3` groups skip this helper because Step 4 assembly already writes directly into `batch_3x3_*`.

## Result Write-Back

Added `write_3x3_batch_solution_to_workspace(...)`, mirroring `write_2x2_batch_solution_to_workspace(...)`.

It copies:

```text
batch_3x3_intercept0_ -> workspace_intercept_[unknown_begin]
batch_3x3_intercept1_ -> workspace_intercept_[unknown_begin + 1]
batch_3x3_intercept2_ -> workspace_intercept_[unknown_begin + 2]
batch_3x3_slope0_     -> workspace_slope_[unknown_begin]
batch_3x3_slope1_     -> workspace_slope_[unknown_begin + 1]
batch_3x3_slope2_     -> workspace_slope_[unknown_begin + 2]
```

This preserves downstream behavior for subtree relation writes and top-down recovery. Those paths are intentionally left for later steps.

## Dense Fallback

The `3x3` dense fallback remains SoA-based. It builds a stack-local dense system from the SoA lane values, solves it with `solve_dense_system(...)`, and writes the result back into the SoA output arrays.

The fallback therefore works for both direct structured SoA assembly and sparse/workspace assembly.

## Behavior Preserved

- Numerical solve formulas are unchanged.
- Sparse-Jacobian validation remains available through `pack_3x3_group_from_workspace(...)`.
- Non-SIMD builds still assemble through workspace and pack into SoA before solving.
- `workspace_intercept_` and `workspace_slope_` are still populated after every `3x3` batch solve.
- Top-down recovery and subtree relation writers were not changed.
- The distributed tree solver was not changed.

## Verification

Commands run after implementation:

```text
cmake --build build/debug --target unittests_reduced_lung --parallel 4
git diff --check
ctest -R "^unittests_reduced_lung$" --output-on-failure
```
