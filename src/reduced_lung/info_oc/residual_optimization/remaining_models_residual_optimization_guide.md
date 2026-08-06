# Remaining Models Residual Optimization Guide

Date: 2026-08-06

Scope: describe the changes needed to extend the optimized residual-evaluation workflow beyond the already optimized common gen16 paths.

## Current Optimized Coverage

The residual optimization so far provides broad low-level improvements plus two model-specific fast paths.

Already applied broadly:

- cached local DOF spans in residual hot loops;
- cached residual value spans in residual hot loops;
- direct local residual writes instead of per-entry `replace_local_value(...)`;
- residual phase profiling;
- serial direct residual norm with distributed global-norm fallback.

Already specialized:

- airway residual: `RigidWall + LinearResistive + no inertia`;
- terminal-unit residual: `LinearElasticity + KelvinVoigt`;
- terminal-unit residual: `LinearElasticity + KelvinVoigt` with all `eta = 0.0`.

Remaining work is to apply the same pattern to the other supported model combinations:

- airway flow models: `LinearResistive`, `NonLinearResistive`;
- airway wall models: `RigidWall`, `KelvinVoigtWall`;
- terminal-unit elasticity models: `LinearElasticity`, `OgdenHyperelasticity`;
- terminal-unit rheology models: `KelvinVoigt`, `FourElementMaxwell`;
- junction and boundary-condition residuals as model-independent residual phases.

## Shared Optimization Pattern

For every remaining model path, use the same workflow:

1. Select fast paths at evaluator construction time, not inside the residual loop.
2. Precompute values that are constant for the whole simulation.
3. Compute dynamic single-use values directly inside the residual consumer loop.
4. Keep persistent model-owned vectors only when values are reused by another phase or by output.
5. Avoid residual-call-time `resize`, clear, fill, copy, and temporary vector readback.
6. Keep generic fallback evaluators for unsupported or less common combinations.
7. Add tests that compare optimized and generic residual values before adding SIMD.

## Airway Residual Extensions

Relevant files:

```text
src/reduced_lung/src/airways/4C_reduced_lung_airways.cpp
src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.cpp
src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.hpp
src/reduced_lung/src/airways/4C_reduced_lung_airways_flow_resistance.cpp
src/reduced_lung/src/airways/4C_reduced_lung_airways_flow_resistance.hpp
```

### RigidWall + LinearResistive + Inertia

Current optimized coverage handles only the no-inertia case.

Needed changes:

- Add a specialized evaluator for `RigidWall + LinearResistive` when any `has_inertia` entry is true.
- Precompute constant Poiseuille resistance once at evaluator construction.
- Precompute constant inertia once at evaluator construction for rigid walls because area is `ref_area`.
- Store zero inertia for elements with `has_inertia == false` so the residual loop stays branch-free.
- Evaluate the complete residual in one direct-write loop.
- Avoid resistance and inertia scratch vectors during residual evaluation.

Expected residual form:

```text
F = p1 - p2 - R * q1 - I / dt * (q1 - q1_n)
```

Implementation note:

- Since `dt` is constant in the current dynamics configuration, optionally precompute `I / dt` during evaluator construction. If future variable time stepping is introduced, keep `I` precomputed and divide by the current `dt` in the loop.

### RigidWall + NonLinearResistive

Needed changes:

- Add specialized evaluators for no-inertia and inertia-enabled variants.
- Precompute the constant Poiseuille base resistance once:

```text
R_poiseuille = 8 * pi * mu * ref_length / ref_area^2
```

- Reuse `NonLinearResistive::k_turb`, which is already model-owned dynamic state.
- Compute the dynamic resistance directly in the residual loop:

```text
R = R_poiseuille * k_turb
```

- If inertia is enabled, reuse the precomputed rigid-wall inertia described above.
- Do not allocate or fill a residual scratch resistance vector.

Important distinction:

- `k_turb` update belongs to the nonlinear-iteration state synchronization path, not the residual callback itself. Optimizing that updater can reduce `state_sync_s`, but it should be handled separately from residual evaluation.

### KelvinVoigtWall + LinearResistive

Current generic path computes resistance and inertia into scratch vectors, then consumes them in the wall residual loop.

Needed changes:

- Add a fused Kelvin-Voigt-wall residual evaluator for `LinearResistive`.
- Compute Poiseuille resistance directly from the current wall area inside the residual loop.
- If inertia is disabled for all elements, skip inertia computation entirely.
- If inertia is enabled, compute the dynamic inertia directly from current area inside the same residual loop.
- Keep using model-owned `viscous_resistance_Rvisc` and `compliance_C`, because those are wall state values reused outside residual evaluation.
- Avoid residual scratch vectors for resistance and inertia.

Expected momentum contribution:

```text
R(area) = 8 * pi * mu * ref_length / area^2
I(area) = rho * ref_length / area
```

### KelvinVoigtWall + NonLinearResistive

Needed changes:

- Add a fused Kelvin-Voigt-wall residual evaluator for `NonLinearResistive`.
- Compute current-area Poiseuille resistance directly in the residual loop.
- Use model-owned `k_turb` directly.
- Compute the nonlinear correction term directly in the same loop.
- If inertia is disabled, avoid all inertia work.
- If inertia is enabled, compute inertia from current area in the same loop.
- Avoid resistance and inertia scratch vectors.

This is the most complex airway residual path. Implement it after validating the simpler Kelvin-Voigt-wall linear-resistance path.

## Terminal-Unit Residual Extensions

Relevant files:

```text
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit.cpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_elasticity.cpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_elasticity.hpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.hpp
```

### LinearElasticity + FourElementMaxwell

Needed changes:

- Add a fused residual evaluator for this pair.
- Compute linear elastic pressure and Four-element Maxwell rheology residual in one loop.
- Update `LinearElasticity::elastic_pressure_p_el` for output compatibility.
- Avoid computing elastic pressure into a temporary vector and then reading it back.
- Precompute constant per-element coefficients involving `dt`, `E_m`, `eta_m`, `eta`, and `reference_volume_v0` if `dt` remains fixed.
- Keep using model-owned `maxwell_pressure_p_m`, because it is a history variable.

Useful precomputed quantities:

```text
denom = E_m * dt + eta_m
branch_viscosity = E_m * dt * eta_m / denom
flow_coeff = (eta + branch_viscosity) / V0
history_coeff = eta_m / denom
```

If future variable time stepping is introduced, recompute the `dt`-dependent coefficients only when `dt` changes.

### OgdenHyperelasticity + KelvinVoigt

Needed changes:

- Add a fused residual evaluator for this pair.
- Compute Ogden elastic pressure and Kelvin-Voigt residual in one loop.
- Update `OgdenHyperelasticity::elastic_pressure_p_el` for output compatibility.
- Add a zero-viscosity variant if all Kelvin-Voigt `eta` entries are `0.0`.
- Avoid elastic-pressure temporary vector readback.

Important caveat:

- Ogden pressure uses `std::pow(...)`, which may dominate this path. First remove temporary-vector overhead, then profile before trying approximate math or SIMD.

### OgdenHyperelasticity + FourElementMaxwell

Needed changes:

- Add a fused residual evaluator for this pair.
- Compute Ogden pressure and Four-element Maxwell residual in one loop.
- Update `OgdenHyperelasticity::elastic_pressure_p_el` for output compatibility.
- Precompute Four-element Maxwell coefficients as described for `LinearElasticity + FourElementMaxwell` when `dt` is fixed.
- Keep model-owned `maxwell_pressure_p_m` as the dynamic history variable.

This path combines the most expensive elasticity model with the more complex rheology model, so profile before and after each small change.

## Junction Residual Extensions

Relevant file:

```text
src/reduced_lung/src/4C_reduced_lung_junctions.cpp
```

Step 6 profiling showed `residual_junctions_s` as the largest residual callback phase after scalar residual optimization.

Needed changes:

- Keep current direct local writes.
- Flatten connection and bifurcation local DOF IDs into SoA arrays during junction construction or local ID assignment.
- Avoid per-junction access through nested containers in the hot residual loop.
- Split connection and bifurcation residual loops into shape-specific kernels.
- Add optional batch/SIMD kernels after scalar SoA flattening is validated.

Suggested SoA fields:

```text
connections_p_out_parent_lid
connections_p_in_child_lid
connections_q_out_parent_lid
connections_q_in_child_lid
connections_first_row

bifurcations_p_out_parent_lid
bifurcations_p_in_child_1_lid
bifurcations_p_in_child_2_lid
bifurcations_q_out_parent_lid
bifurcations_q_in_child_1_lid
bifurcations_q_in_child_2_lid
bifurcations_first_row
```

This is a good SIMD candidate because the residual equations are simple and repeated many times.

## Boundary-Condition Residual Extensions

Relevant file:

```text
src/reduced_lung/src/4C_reduced_lung_boundary_conditions.cpp
```

Current behavior already evaluates function-valued BCs once per BC model and writes directly.

Needed changes:

- Keep direct local writes.
- Split constant-value and function-value BC loops into simple shape-specific kernels.
- For constant pressure or flow BCs, ensure values are stored in contiguous arrays and avoid per-entry indirection where possible.
- For function BCs, keep evaluating the function once per model per residual evaluation.
- Add batch/SIMD only if `residual_boundary_conditions_s` remains significant after junction and model-specific residual work.

## SIMD and Batch Strategy

Do not add SIMD until scalar specialized paths have been implemented and reprofiled.

When SIMD is justified, follow this order:

1. Junction residuals, because Step 6 showed they are the largest residual callback phase.
2. Rigid airway residuals, because their formulas are simple and data is already SoA-like.
3. Linear terminal-unit residuals, especially `LinearElasticity + KelvinVoigt` and `LinearElasticity + FourElementMaxwell`.
4. Boundary-condition residuals if still visible.
5. Ogden paths only after measuring whether `std::pow(...)` dominates.

SIMD requirements:

- group by model shape and equation shape;
- use SoA arrays;
- use masked or padded tails;
- keep scalar fallback paths;
- validate residual vectors against scalar reference paths.

## Validation Requirements

For each new fast path, add a direct-vs-generic residual comparison test before relying on performance numbers.

Suggested tests:

```text
RigidLinearInertiaResidualMatchesGenericPath
RigidNonlinearResidualMatchesGenericPath
KelvinVoigtWallLinearResidualMatchesGenericPath
KelvinVoigtWallNonlinearResidualMatchesGenericPath
LinearFourElementMaxwellResidualMatchesGenericPath
OgdenKelvinVoigtResidualMatchesGenericPath
OgdenFourElementMaxwellResidualMatchesGenericPath
JunctionSoAResidualMatchesGenericPath
BoundaryConditionBatchResidualMatchesGenericPath
```

Run focused checks after each path:

```bash
cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_(tree|sparse).*\.4C\.yaml-p1$" --output-on-failure
```

For distributed-sensitive changes, also run:

```bash
cmake --build build/debug --target unittests_reduced_lung.np2 --parallel 4
ctest -R "^unittests_reduced_lung\.np2$" --output-on-failure
```

## Recommended Implementation Order

1. Add `RigidWall + LinearResistive + inertia` fast path.
2. Add `RigidWall + NonLinearResistive` fast paths.
3. Add `LinearElasticity + FourElementMaxwell` fused terminal-unit path.
4. Add `OgdenHyperelasticity + KelvinVoigt` fused terminal-unit path.
5. Add `OgdenHyperelasticity + FourElementMaxwell` fused terminal-unit path.
6. Add `KelvinVoigtWall + LinearResistive` fused airway path.
7. Add `KelvinVoigtWall + NonLinearResistive` fused airway path.
8. Flatten junction residual metadata to SoA and reprofile.
9. Consider SIMD/batch kernels only for phases still significant after scalar specialization.

## Expected Outcome

The optimized residual workflow should eventually have these properties for every supported model combination:

- evaluator construction chooses the fastest valid residual path;
- constant coefficients are precomputed once;
- dynamic single-use values are computed in the consumer loop;
- persistent vectors are updated only when another phase or output needs them;
- residual loops use cached spans and direct local writes;
- generic fallback paths remain correct;
- SIMD is added only where profiling still justifies it.
