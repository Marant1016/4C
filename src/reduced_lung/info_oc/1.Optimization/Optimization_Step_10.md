# Optimization Step 10 - Direct Tree-Block Coefficient Storage

## Scope

This step implements the direct coefficient storage follow-up after Step 9.

Only the serial `TreeNewtonLinearSolver` coefficient-access path was changed. The distributed `DistributedTreeNewtonLinearSolver` implementation was intentionally left unchanged.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Added Direct Coefficient Value Buffers

The serial solver now owns numeric coefficient buffers that mirror the Step 9 coefficient-location descriptor arrays:

```cpp
double root_boundary_coefficient_value_ = 0.0;
std::vector<double> equation_inlet_pressure_coefficient_values_;
std::vector<double> matrix_coefficient_values_;
std::vector<double> child_pressure_parent_coefficient_values_;
std::vector<double> child_pressure_child_coefficient_values_;
std::vector<double> child_flow_coefficient_values_;
```

The dense local matrix coefficient values use the same indexing as `workspace_matrix_`. The inlet-pressure coefficient values use the same per-equation indexing as `workspace_rhs_inlet_pressure_`.

### Populate Direct Values Per Structured Linearization

`resolve_structured_coefficient_locations(...)` now also copies coefficient values from the current `TreeLinearization` into the direct buffers.

Missing optional coefficients are copied as `0.0`, preserving the previous `TreeLinearization::value(...)` behavior.

### Use Direct Values During Structured Solves

The serial structured coefficient provider now has an overload that accepts both the symbolic descriptor and the direct value.

For `StructuredTreeBlocks`, solve-time bottom-up assembly reads the direct value and no longer dereferences `TreeLinearization` row entries.

For `SparseJacobian`, the same call path ignores the direct value and uses the descriptor's `(local_row, local_dof)` pair to perform the existing sparse row scan. This keeps sparse validation behavior unchanged.

Required coefficients still pass through `required_matrix_value(...)`, so missing or near-zero required coefficients still fail loudly.

## Behavior Preserved

- The bottom-up condensation formula is unchanged.
- The top-down recovery formula is unchanged.
- Missing optional structured coefficients still behave as zero.
- Required coefficient assertions are unchanged.
- Sparse-Jacobian coefficient access still uses the existing sparse row scan.
- The distributed tree solver was not changed.

## Intended Benefit

Step 9 moved repeated structured row-entry searches out of the solve loop. Step 10 goes one step further: the structured solve loop now reads compact coefficient arrays directly instead of reading through `TreeLinearization::entries(...)`.

This should reduce solve-time `tree_lookup_s` for the serial `StructuredTreeBlocks` path. The copy cost is paid when `set_tree_linearization(...)` receives a fresh structured linearization, which is part of structured tree linearization assembly in the custom Newton path.

## Verification

Verification passed:

```text
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

The reduced-lung unit test and focused `newton_tree` input tests passed.

No release benchmark rerun was performed in this step.

## Follow-Up

The next release benchmark run should compare Step 10 against Step 9 and Step 8. The important counters are `tree_lookup_s`, `tree_solve_s`, `tree_assembly_s`, `linear_solve_s`, and full-solve `newton_total_s`, because this step intentionally shifts coefficient access work away from solve-time row-entry reads.
