# Optimization Step 14 - Reuse Structured Linearization And Coefficient Locations

## Scope

This step reduces repeated structured-tree setup work during the serial `NewtonTree` workflow.

Only the custom Newton structured-tree assembly path and the serial `TreeNewtonLinearSolver` coefficient-location refresh were changed. The distributed tree solver was intentionally left unchanged.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_newton_solver.cpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Reused TreeLinearization Storage

`NewtonSolver::assemble_tree_linearization_for_current_state()` now reuses the existing `TreeLinearization` storage when the row and dof dimensions are unchanged:

```cpp
if (tree_linearization_.num_rows() == num_rows && tree_linearization_.num_dofs() == num_dofs)
{
  tree_linearization_.clear_values();
}
else
{
  tree_linearization_.reset(num_rows, num_dofs);
}
```

This keeps per-row vector capacity between Newton corrections instead of rebuilding the row storage every time.

### Reused Structured Coefficient Entry Indices

`TreeNewtonLinearSolver::resolve_structured_coefficient_locations()` now first tries the previously resolved `structured_entry_index`.

If the cached index is valid and still points to the expected dof, the solver reuses it directly. If not, it falls back to the existing row scan and updates the cached index.

## Behavior Preserved

- The Newton equations are unchanged.
- The tree condensation and top-down recovery mathematics are unchanged.
- Missing coefficients still resolve to `0.0` as before.
- Existing row and dof assertions are preserved.
- Sparse-Jacobian solver behavior is unchanged.
- The distributed tree solver was not changed.

## Expected Performance Effect

This step targets full-solve `tree_assembly_s` by reducing repeated allocation and repeated coefficient-location scans between Newton corrections.

The biggest expected benefit is for full `NewtonTree` solves where structured tree assembly is now comparable to, or larger than, the tree linear solve itself.

## Verification

Verification passed:

```text
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
git diff --check
```

The `ctest` commands were run from `build/debug`.

Results:

- `unittests_reduced_lung`: passed, `1/1` tests.
- Focused `newton_tree` input tests: passed, `7/7` tests.
- `git diff --check`: passed.

## Benchmark Status

Release benchmarks were not run for this step.

When benchmarking later, compare against Step 13 with focus on:

- `tree_assembly_s`
- `newton_total_s`
- `tree_solve_s`
- nonlinear iteration count

## Conclusion

Step 14 keeps the same solver behavior but avoids repeated structured-tree setup work. It should be kept if future release benchmarks show lower full-solve tree assembly time without increasing tree solve time or nonlinear iterations.
