# Implementation NextStep 3 - Reduce Allocations And Reuse Symbolic Information

This step keeps the existing serial tree elimination algorithm and structured coefficient source, but separates topology-dependent setup from per-Newton-solve numeric updates.

Runtime `NewtonTree` remains serial-only and still consumes structured tree-linearization coefficients assembled by the reduced-lung physics modules.

## Tree Solver Refactor

Updated:

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

`TreeNewtonLinearSolver` now builds a symbolic solve plan in its constructor. The plan caches per-element data that does not change between Newton corrections:

- unknown global dof ids
- unknown local dof ids
- local equation row ids
- inlet/outlet dof lookup indices
- child-interface rows and dof ids
- parent-outlet pressure recovery index
- root inlet boundary row and dof id

The per-solve path now reuses solver-owned work buffers:

- one dense local matrix per tree element
- two local right-hand sides per tree element
- per-element intercept and inlet-pressure slope vectors
- bottom-up subtree relations
- top-down inlet-pressure corrections
- per-interface child pressure slope/intercept values

## Numeric Solve Changes

The old local solve path assembled a fresh dense matrix/RHS set for each element and solved two systems by copying the matrix twice:

```text
local_matrix * intercept = constant_rhs
local_matrix * slope     = inlet_pressure_rhs
```

The new path refills the preallocated local matrix and both RHS vectors, then performs one Gaussian elimination for both RHS vectors. The dense matrix and RHS buffers are intentionally overwritten and reused on the next solve.

## Interface-Recovery Reuse

During bottom-up condensation, the solver already computes each child inlet pressure relation:

```text
child_inlet_pressure = slope * parent_outlet_pressure + intercept
```

Those per-interface numeric relations are now stored in the element workspace. The top-down recovery pass reuses them instead of looking up the same pressure-continuity coefficients again.

## Test Updates

Extended:

- `src/reduced_lung/tests/4C_reduced_lung_tree_linear_solver_test.cpp`

Added `ReusedStructuredTreeSolverMatchesSparseSolver`, which reuses one structured-source `TreeNewtonLinearSolver` instance across several mixed airway/terminal-unit state and time points. Each iteration assembles a fresh sparse reference system and structured tree linearization, then checks that the reused tree solver correction matches the sparse linear-solver correction.

## Verification Run

Commands run successfully:

```text
cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_.*\.4C\.yaml-p1$" --output-on-failure
cmake --build build/debug --target unittests_reduced_lung.np2 --parallel 4
ctest -R "^unittests_reduced_lung(\.np2)?$" --output-on-failure
git diff --check
```

The runtime YAML regex also selected `test_cleanup`; all selected tests passed.

## Remaining After This Step

- `NewtonTree` is still serial-only.
- Sparse Jacobian support remains required for `Nox`, `NewtonSparse`, and sparse-reference tree tests.
- Runtime `NewtonTree` still assembles structured numeric coefficients every Newton correction.
- Performance benchmarking/profiling is still future work.
- Parallel tree solving is still future work.
