# Implementation NextStep 1 - Validate The Existing Tree Path

This step adds validation for the already implemented `NewtonTree` workflow. It does not replace the current sparse-Jacobian coefficient source and does not make `NewtonTree` the default.

## Code Changes

- Extended `src/reduced_lung/tests/4C_reduced_lung_tree_linear_solver_test.cpp`.
- Added sparse-vs-tree Newton correction checks for additional model combinations:
  - Kelvin-Voigt airway wall model.
  - Nonlinear airway resistance.
  - Four-element Maxwell terminal unit rheology.
  - Larger mixed tree with airways and terminal units.
- Added serial full time-step workflow comparisons using the same in-memory problem setup for:
  - `Nox`.
  - `NewtonSparse`.
  - `NewtonTree`.
- The full workflow comparisons check solution vectors, owned dofs, residual norms, selected flow balances, and terminal-unit volume updates.
- Added optional nonzero initial-state seeding in the tree linear-solver test fixture so nonlinear resistance and Kelvin-Voigt airway derivatives are exercised away from the trivial zero-flow state.

## Runtime YAML Validation

Added serial runtime input variants for representative reduced-lung cases:

- `tests/input_files/reduced_lung_terminal_unit_newton_sparse.4C.yaml`
- `tests/input_files/reduced_lung_terminal_unit_newton_tree.4C.yaml`
- `tests/input_files/reduced_lung_serial_airways_flow_newton_sparse.4C.yaml`
- `tests/input_files/reduced_lung_serial_airways_flow_newton_tree.4C.yaml`
- `tests/input_files/reduced_lung_aw_bifurcation_flow_newton_sparse.4C.yaml`
- `tests/input_files/reduced_lung_aw_bifurcation_flow_newton_tree.4C.yaml`
- `tests/input_files/reduced_lung_3_aw_2_tu_newton_sparse.4C.yaml`
- `tests/input_files/reduced_lung_3_aw_2_tu_newton_tree.4C.yaml`
- `tests/input_files/reduced_lung_3_aw_2_tu_4elemax_and_kv_newton_sparse.4C.yaml`
- `tests/input_files/reduced_lung_3_aw_2_tu_4elemax_and_kv_newton_tree.4C.yaml`

Registered those files in `tests/list_of_tests.cmake` as serial `p1` CTests. The original YAML files still omit `nonlinear_solver`, so they continue to use the default `Nox` workflow.

## Supporting Fix

The new Kelvin-Voigt airway validation exposed that the reduced-lung sparse Jacobian allocation used an estimate of 3 entries per row, while Kelvin-Voigt airway state rows insert 4 coefficients.

Updated sparse Jacobian construction to use 4 estimated entries per row in:

- `src/reduced_lung/src/4C_reduced_lung_main.cpp`
- the reduced-lung tree linear-solver test fixture

This is not a tree-solver algorithm change; it allows existing Kelvin-Voigt airway Jacobian assembly to run through the reduced-lung sparse matrix path.

## Solver Policy

Explicit `NewtonTree` remains strict. Unsupported tree layouts or unsupported MPI usage still throw. No silent fallback to `NewtonSparse` was added.

## Verification Run

Commands run successfully:

```text
cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_.*\.4C\.yaml-p1$" --output-on-failure
cmake --build build/debug --target unittests_reduced_lung.np2 --parallel 4
ctest -R "^unittests_reduced_lung(\.np2)?$" --output-on-failure
```

The new runtime YAML test command also matched and ran `test_cleanup` as part of the regex selection; all selected tests passed.

## Remaining After This Step

- `NewtonTree` is still serial-only.
- `NewtonTree` still extracts coefficients from the assembled sparse Jacobian.
- Sparse-Jacobian-free structured derivative/block assembly remains future work.
- Larger production-scale tree benchmarks are still not implemented.
