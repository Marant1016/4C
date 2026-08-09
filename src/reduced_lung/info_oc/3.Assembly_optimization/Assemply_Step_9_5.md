# Assemply Step 9.5

## Goal

Optimize the Ogden terminal-unit elasticity pressure-gradient loop without changing the model equation.

This step touches only terminal-unit elasticity gradient evaluation. It does not change rheology, tree solver logic, coefficient write APIs, or sparse assembly behavior.

## Target

File:

```text
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_elasticity.cpp
```

Function:

```text
ogden_hyperelastic_pressure_gradient(...)
```

## What Changed

Before this step, the Ogden gradient loop repeatedly accessed shared model and data vectors inside every element iteration:

```text
locally_relevant_dofs.local_values_as_span()
data.lid_q
data.volume_v
data.reference_volume_v0
ogden_hyperelastic_model.bulk_modulus_kappa
ogden_hyperelastic_model.nonlinear_stiffening_beta
ogden_hyperelastic_model.elastic_pressure_grad_dp_el
```

After this step, the loop now:

```text
1. Stores the local dof span once before the loop.
2. Hoists references to terminal-unit data vectors and Ogden model vectors.
3. Stores per-element reference volume and beta values locally.
4. Computes v0_over_vi once per element.
5. Reuses v0_over_vi^2 and pow(v0_over_vi, beta) in the final expression.
```

The pressure-gradient formula is unchanged:

```text
dp_el/dq = kappa * dt / (beta * v0)
           * (v0 / vi)^2
           * ((beta + 1) * (v0 / vi)^beta - 1)
```

where:

```text
vi = volume + dt * q
```

## Why

Ogden terminal units require dynamic pressure-gradient values during nonlinear assembly. The dominant `std::pow(...)` call remains, but the surrounding repeated vector/span access and duplicated intermediate calculations are reduced.

## Expected Effect

Expected improvement is limited if `std::pow(...)` dominates the loop, but this still removes avoidable overhead for inputs with many Ogden terminal units.

The gen16 rigid/linear baseline is not expected to show a meaningful timing change from this step unless Ogden terminal units are active.

## Verification

Built the reduced-lung unit target:

```bash
cmake --build /scratch/Rodriguez/workspace/4C/4C/build/release \
  --target unittests_reduced_lung \
  --parallel 4
```

Ran reduced-lung unit tests:

```bash
ctest -R "^unittests_reduced_lung$" --output-on-failure
```

Result:

```text
100% tests passed, 0 tests failed out of 1
```

## Safety

Step 9.1 added direct-vs-generic structured assembly coverage for Ogden terminal-unit paths. That coverage protects this optimized gradient loop by checking that direct structured coefficients still match the generic tree linearization path and that the direct tree correction still matches the sparse solver.
