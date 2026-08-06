# Residual Step 10: Junction and Boundary-Condition Residual Phases

Date: 2026-08-06

## Scope

Implemented the remaining non-model residual-phase optimizations from `remaining_models_residual_optimization_guide.md`:

- junction residual SoA metadata and split residual kernels;
- boundary-condition residual batching for constant and function-valued models.

## Junction Changes

Updated:

```text
src/reduced_lung/src/4C_reduced_lung_junctions.hpp
src/reduced_lung/src/4C_reduced_lung_junctions.cpp
```

Changes:

- Added flattened local row and local DOF arrays for connection residuals.
- Added flattened local row and local DOF arrays for bifurcation residuals.
- Populated the flattened arrays during local equation and local DOF assignment.
- Rewrote junction residual assembly to use split connection and bifurcation loops over the flattened arrays.
- Kept the existing nested local DOF arrays for Jacobian, tree-linearization, and compatibility with existing setup code.

## Boundary-Condition Changes

Updated:

```text
src/reduced_lung/src/4C_reduced_lung_boundary_conditions.cpp
```

Changes:

- Changed boundary-condition residual assembly to batch constant-value models first.
- Changed boundary-condition residual assembly to batch function-valued models separately.
- Kept function-valued BC behavior of evaluating the function once per model per residual evaluation.
- Kept direct local residual writes.

## Validation Tests Added

Updated:

```text
src/reduced_lung/tests/4C_reduced_lung_junctions_test.cpp
src/reduced_lung/tests/4C_reduced_lung_boundary_conditions_test.cpp
```

Added focused checks:

- `JunctionSoAResidualMatchesGenericPath`
- `BoundaryConditionBatchResidualMatchesGenericPath`

Expanded local DOF assignment checks to verify the new flattened junction metadata.

## Verification

Executed focused checks only:

```bash
cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_(tree|sparse).*\.4C\.yaml-p1$" --output-on-failure
cmake --build build/debug --target unittests_reduced_lung.np2 --parallel 4
ctest -R "^unittests_reduced_lung\.np2$" --output-on-failure
```

Results:

- `unittests_reduced_lung` build passed.
- `unittests_reduced_lung` passed.
- Reduced-lung Newton tree/sparse input tests passed: 13/13.
- `unittests_reduced_lung.np2` build passed.
- `unittests_reduced_lung.np2` passed.

## Unresolved Issues

None for the residual optimization items described in the guide.
