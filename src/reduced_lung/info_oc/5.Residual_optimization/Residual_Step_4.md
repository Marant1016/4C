# Residual Optimization Step 4

Date: 2026-08-06

## Scope

Implemented Step 4 from `residual_evaluation_optimization_guide.md`: fuse the common terminal-unit residual for `Linear elasticity + Kelvin-Voigt rheology`. No other terminal-unit residual combinations were changed.

## Implemented Changes

- Added a fused residual evaluator for `LinearElasticity` with `KelvinVoigt` rheology.
- Added a dedicated zero-viscosity variant for the gen16-style `eta = 0` case.
- The fused evaluator computes elastic pressure and residual in one loop.
- The residual uses the local elastic-pressure value directly instead of reading it back from the stored pressure vector.
- The stored `elastic_pressure_p_el` vector is still updated so high-verbosity output keeps the same data source.
- Existing generic residual evaluation remains active for:
  - Ogden elasticity with Kelvin-Voigt rheology;
  - Four-element Maxwell rheology;
  - future unsupported combinations.

## Modified Files

- `src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit.cpp`
- `src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp`
- `src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.hpp`

## Diff Summary

Focused source diff:

```text
3 files changed, 89 insertions(+), 6 deletions(-)
```

Change distribution:

- `4C_reduced_lung_terminal_unit.cpp`: passes the elasticity model into the rheology residual evaluator factory.
- `4C_reduced_lung_terminal_unit_rheology.hpp`: updates the residual evaluator factory signature.
- `4C_reduced_lung_terminal_unit_rheology.cpp`: adds fused linear Kelvin-Voigt residual evaluators and selects them when applicable.

## Verification

Build:

```bash
cmake --build build/debug --target unittests_reduced_lung --parallel 4
```

Result: passed.

Tests:

```bash
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_(terminal_unit|3_aw_2_tu).*\.4C\.yaml-p1$" --output-on-failure
```

Result: both passed.

## Notes

- The zero-viscosity branch is selected only when every Kelvin-Voigt viscosity entry is exactly `0.0`.
- The generic pressure-evaluator path is retained for non-linear elasticity and Four-element Maxwell cases.
- No Jacobian, tree-linearization, state-update, output, or SIMD code was changed in this step.
