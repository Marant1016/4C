# Residual Optimization Step 3

Date: 2026-08-06

## Scope

Implemented Step 3 from `residual_evaluation_optimization_guide.md`: specialize the common airway residual for the `Rigid wall + Linear resistance + no inertia` model combination. No other residual algorithms were changed.

## Implemented Changes

- Added a direct residual evaluator for rigid-wall, linear-resistance airways without inertia.
- Added a helper to detect whether any airway element in a `LinearResistive` model has inertia enabled.
- Added lazy one-time precomputation of constant Poiseuille resistance for the specialized path.
- The specialized evaluator computes the full residual in one loop and writes directly into local residual storage.
- The specialized path avoids the previous per-evaluation resistance scratch fill and inertia scratch fill.
- Generic residual evaluation remains active for:
  - nonlinear resistance;
  - any inertia-enabled linear resistance model;
  - Kelvin-Voigt wall models;
  - unsupported or future model combinations.

## Modified Files

- `src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.cpp`

## Diff Summary

Focused source diff:

```text
1 file changed, 61 insertions(+)
```

The diff adds:

- `<numbers>` for the Poiseuille factor using `std::numbers::pi`;
- inertia-detection and resistance-precompute helpers;
- the direct rigid-linear-no-inertia residual evaluator;
- evaluator selection logic in `make_residual_evaluator(...)`.

## Verification

Build:

```bash
cmake --build build/debug --target unittests_reduced_lung --parallel 4
```

Result: passed.

Tests:

```bash
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_(serial_airways_flow|aw_bifurcation_flow|lung_tree_gen10|3_aw_2_tu).*\.4C\.yaml-p1$" --output-on-failure
```

Result: both passed.

## Notes

- The specialization is selected only when the wall model is `RigidWall`, the flow model is `LinearResistive`, and all `has_inertia` entries are false.
- The precomputed resistance vector is retained by the residual-evaluator lambda and reused across residual evaluations.
- No SIMD path was added in this step.
- No terminal-unit, junction, boundary-condition, tree-linearization, or solver logic was modified in this step.
