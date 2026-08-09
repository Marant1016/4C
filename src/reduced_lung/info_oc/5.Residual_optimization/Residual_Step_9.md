# Residual Step 9: Remaining Model Fast Paths

Date: 2026-08-06

## Scope

Implemented the remaining scalar model-combination residual fast paths described in `remaining_models_residual_optimization_guide.md`.

This step covers model-specific residual paths only. Junction SoA flattening and boundary-condition batching are model-independent residual-phase work and were not changed here.

## Implemented Airway Changes

Updated:

```text
src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.cpp
```

Changes:

- Added `RigidWall + NonLinearResistive` fast paths.
- Reused precomputed rigid Poiseuille resistance as the nonlinear base resistance.
- Reused precomputed rigid inertia for inertia-enabled nonlinear rigid airways.
- Added fused `KelvinVoigtWall + LinearResistive` residual evaluation.
- Added fused `KelvinVoigtWall + NonLinearResistive` residual evaluation.
- Removed residual-call scratch resistance and inertia buffers from the optimized airway model paths.
- Kept generic fallback evaluators for unsupported future combinations.

## Implemented Terminal-Unit Changes

Updated:

```text
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit.cpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.hpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp
```

Changes:

- Added fused `LinearElasticity + FourElementMaxwell` residual evaluation.
- Added fused `OgdenHyperelasticity + KelvinVoigt` residual evaluation.
- Added zero-viscosity `OgdenHyperelasticity + KelvinVoigt` residual evaluation.
- Added fused `OgdenHyperelasticity + FourElementMaxwell` residual evaluation.
- Updated elastic-pressure output state inside each fused path.
- Added construction-time allocation of Four-element Maxwell residual coefficient buffers.
- Cached Four-element Maxwell `flow_coeff` and `history_coeff` values per time-step size and recomputed them only if `dt` changes.
- Kept generic pressure-evaluator fallback paths.

## Validation Tests Added

Updated:

```text
src/reduced_lung/tests/4C_reduced_lung_airways_test.cpp
src/reduced_lung/tests/4C_reduced_lung_terminal_unit_test.cpp
```

Added focused residual checks:

- `RigidNonlinearResidualMatchesGenericPath`
- `KelvinVoigtWallLinearResidualMatchesGenericPath`
- `KelvinVoigtWallNonlinearResidualMatchesGenericPath`
- `LinearFourElementMaxwellResidualMatchesGenericPath`
- `OgdenKelvinVoigtResidualMatchesGenericPath`
- `OgdenFourElementMaxwellResidualMatchesGenericPath`

The existing Step 8 test continues to cover `RigidWall + LinearResistive + inertia`.

## Verification

Executed focused checks only:

```bash
cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_(tree|sparse).*\.4C\.yaml-p1$" --output-on-failure
```

Results:

- `unittests_reduced_lung` build passed.
- `unittests_reduced_lung` passed.
- Reduced-lung Newton tree/sparse input tests passed: 13/13.

## Unresolved Issues

- Junction SoA flattening remains separate model-independent residual-phase work.
- Boundary-condition batching remains separate model-independent residual-phase work.
