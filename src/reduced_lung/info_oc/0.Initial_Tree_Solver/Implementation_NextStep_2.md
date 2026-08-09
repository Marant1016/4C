# Implementation NextStep 2 - Move Toward Structured Local Blocks

This step moves the runtime `NewtonTree` path away from generic sparse-Jacobian coefficient extraction. The tree solver can now consume structured local derivative blocks assembled directly by the reduced-lung physics modules.

The old sparse coefficient source is still available as a reference path for tests.

## New Structured Linearization Storage

Added:

- `src/reduced_lung/src/4C_reduced_lung_tree_linearization.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linearization.cpp`

Main type:

```text
TreeLinearization
```

`TreeLinearization` stores coefficients addressed by:

```text
local residual row id
local locally-relevant dof id
```

This mirrors the coefficient lookup previously done through `SparseMatrix::extract_my_row_view(...)`, but without requiring a completed global sparse Jacobian.

## Assembly Pipeline Extension

Extended `ReducedLungAssemblyPipeline` with:

```text
tree_linearization_assemblers
```

These callbacks are registered next to the existing residual and sparse-Jacobian callbacks in:

```text
src/reduced_lung/src/4C_reduced_lung_helpers.cpp
```

The existing sparse Jacobian callbacks remain unchanged and are still used by `Nox`, `NewtonSparse`, and sparse-reference tree tests.

## Structured Physics Blocks

Added structured derivative/block assembly in the existing physics modules:

- airways
- terminal units
- junctions
- boundary conditions

Updated files include:

- `src/reduced_lung/src/airways/4C_reduced_lung_airways_common.hpp`
- `src/reduced_lung/src/airways/4C_reduced_lung_airways.hpp`
- `src/reduced_lung/src/airways/4C_reduced_lung_airways.cpp`
- `src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.hpp`
- `src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.cpp`
- `src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_common.hpp`
- `src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit.hpp`
- `src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit.cpp`
- `src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.hpp`
- `src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp`
- `src/reduced_lung/src/4C_reduced_lung_junctions.hpp`
- `src/reduced_lung/src/4C_reduced_lung_junctions.cpp`
- `src/reduced_lung/src/4C_reduced_lung_boundary_conditions.hpp`
- `src/reduced_lung/src/4C_reduced_lung_boundary_conditions.cpp`

The derivative formulas remain in their physics modules. The tree solver only consumes the resulting coefficients.

## Tree Solver Refactor

Updated:

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

`TreeNewtonLinearSolver` now supports two coefficient sources:

```text
SparseJacobian
StructuredTreeBlocks
```

The bottom-up/top-down tree elimination algorithm is unchanged. Only coefficient access was abstracted.

Sparse reference behavior:

```text
SparseMatrix -> coefficient lookup -> tree solve
```

Structured runtime behavior:

```text
TreeLinearization -> coefficient lookup -> tree solve
```

## Newton Runtime Change

Updated:

- `src/reduced_lung/src/4C_reduced_lung_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_newton_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_newton_solver.cpp`
- `src/reduced_lung/src/4C_reduced_lung_main.cpp`

`NewtonLinearSolver` can now report which linearization it needs:

```text
SparseJacobian
StructuredTreeBlocks
```

Runtime `NewtonTree` now constructs `TreeNewtonLinearSolver` with `StructuredTreeBlocks`.

The custom Newton loop now does:

```text
Nox:
  unchanged NOX sparse Jacobian path

NewtonSparse:
  residual assembly
  sparse Jacobian assembly
  sparse linear solve

NewtonTree:
  residual assembly
  structured tree-linearization assembly
  tree solve
```

So runtime `NewtonTree` no longer assembles or completes the sparse Jacobian for Newton corrections.

## Test Updates

Extended:

```text
src/reduced_lung/tests/4C_reduced_lung_tree_linear_solver_test.cpp
```

The direct correction tests now compare:

```text
SparseNewtonLinearSolver correction
sparse-source TreeNewtonLinearSolver correction
structured-source TreeNewtonLinearSolver correction
```

The full time-step workflow tests now use structured-source `NewtonTree`, matching runtime behavior.

The NextStep 1 YAML runtime tests for `NewtonTree` also now exercise the structured tree path.

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

The YAML regex also selected `test_cleanup`; all selected tests passed.

## Remaining After This Step

- `NewtonTree` is still serial-only.
- Sparse Jacobian support remains required for `Nox`, `NewtonSparse`, and sparse-reference tree tests.
- Structured tree blocks are now assembled every Newton correction; symbolic/preallocated block plans are still future work.
- Performance benchmarking/profiling is still future work.
- Parallel tree solving is still future work.
