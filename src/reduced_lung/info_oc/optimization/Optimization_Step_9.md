# Optimization Step 9 - Precompute Structured Coefficient Locations

## Scope

This step implements the first follow-up optimization from the 2026-07-31 serial tree solver guide.

Only the serial `TreeNewtonLinearSolver` coefficient-access path was changed. The distributed `DistributedTreeNewtonLinearSolver` implementation was intentionally left unchanged.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Added Coefficient Location Descriptors

The serial solver now stores symbolic coefficient requests as compact descriptors:

```cpp
struct TreeCoefficientLocation
{
  int local_row = -1;
  int local_dof = -1;
  int structured_entry_index = -1;
};
```

The descriptors keep the existing local residual row and local dof ids, plus a resolved row-entry index for `TreeLinearization` values.

### Precomputed Serial Solver Coefficient Requests

`build_symbolic_plan()` now creates coefficient-location arrays for the hot serial tree-solve accesses:

- root inlet boundary coefficient,
- inlet-pressure coefficient for each element equation row,
- dense local matrix coefficient for each element block entry,
- pressure-continuity parent-pressure coefficient for each child interface,
- pressure-continuity child-pressure coefficient for each child interface,
- junction child-flow coefficient for each child interface.

The dense matrix coefficient descriptors use the same indexing as `workspace_matrix_`, and the inlet-pressure equation descriptors use the same per-row indexing as `workspace_rhs_inlet_pressure_`.

### Resolved Structured Entry Indices Once Per Linearization

`set_tree_linearization(...)` now resolves each descriptor's `structured_entry_index` when the serial solver is using `StructuredTreeBlocks`.

This moves the row search out of the tree solve hot path. Missing optional coefficients keep `structured_entry_index = -1` and still read as zero.

### Switched Hot Reads To Descriptor-Based Access

The bottom-up assembly and root-boundary setup now read coefficients through the precomputed descriptors.

For the structured provider, descriptor-based reads access:

```cpp
tree_linearization.entries(location.local_row)[location.structured_entry_index].second
```

For the sparse provider, descriptor-based reads fall back to the existing sparse row scan using `location.local_row` and `location.local_dof`. This keeps sparse validation behavior unchanged.

Required coefficients still use `required_matrix_value(...)`, so missing or near-zero required entries still fail loudly.

## Behavior Preserved

- The bottom-up condensation formula is unchanged.
- The top-down recovery formula is unchanged.
- Missing optional matrix coefficients still behave as zero.
- Required coefficients still trigger the existing missing/near-zero assertion.
- Sparse-Jacobian coefficient access still uses the existing sparse row scan.
- The distributed tree solver was not changed.

## Intended Benefit

This step targets `tree_lookup_s`, which remained a significant cost after Step 8. The structured serial solver no longer repeats a linear search through `TreeLinearization` row entries for every hot-path coefficient access during bottom-up assembly.

The lookup work is now paid once per structured tree linearization in `set_tree_linearization(...)`, which is called after the structured coefficients are assembled for the current Newton state.

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

The next benchmark run should compare Step 9 against Step 8 with longer release benchmark timing, focusing on `tree_lookup_s`, `tree_solve_s`, `linear_solve_s`, and full-solve `newton_total_s`.
