# Assemply Step 9.3

## Goal

Optimize Kelvin-Voigt airway resistance and inertia derivative loops without changing the model equations.

This step touches only airway derivative evaluation. It does not change tree solver logic, batching policy, or coefficient write APIs.

## Target

File:

```text
src/reduced_lung/src/airways/4C_reduced_lung_airways_flow_resistance.cpp
```

Functions:

```text
evaluate_linear_flow_resistance_derivative_kelvin_voigt(...)
evaluate_nonlinear_flow_resistance_derivative_kelvin_voigt(...)
evaluate_inertia_derivative_kelvin_voigt(...)
```

## What Changed

The previous loops repeatedly accessed values such as:

```text
dofs.local_values_as_span()
data.air_properties.density
data.air_properties.dynamic_viscosity
area[i] * area[i]
area[i] * area[i] * area[i]
q1 + q2
q2 - q1
```

The updated loops now:

```text
1. Store the local dof span once before each loop.
2. Hoist density, viscosity, pi, and dt factors where possible.
3. Store per-element values such as area, area^2, area^3, ref_length, q1, and q2 locally.
4. Reuse q_sum and q_difference variables in the nonlinear Kelvin-Voigt derivative.
5. Keep turbulence and inertia branch logic unchanged.
```

## Linear Kelvin-Voigt Resistance Derivative

The linear Kelvin-Voigt derivative still computes the same two values:

```text
resistance_derivative_q1
resistance_derivative_q2
```

but now avoids repeated span access and repeated area-power calculations.

## Nonlinear Kelvin-Voigt Resistance Derivative

The nonlinear derivative still uses:

```text
k_turb
dk_dq1 / dk_dq2
alpha
dalpha_dk
q_sum = q1 + q2
q_difference = q2 - q1
```

The formulas are unchanged. The code now computes the reused terms once per element and reuses them in both q1 and q2 derivative expressions.

## Kelvin-Voigt Inertia Derivative

The inertia derivative still respects the existing `has_inertia` branch:

```text
if inertia is disabled:
  inertia_q1 = 0
  inertia_q2 = 0

if inertia is enabled:
  compute q1/q2 inertia derivatives
```

The optimized loop stores the dof span, density, and `-0.5 / dt` factor once, then reuses local per-element values.

## Why

Kelvin-Voigt airway blocks produce 3x3 structured tree coefficients. These derivative loops are part of the dynamic coefficient assembly before the tree solve.

The coefficients still have to be recomputed each Newton iteration, but avoidable overhead inside the loop can be reduced.

## Expected Effect

Expected improvement is small to medium for input files with many Kelvin-Voigt airways.

The exact gain depends on whether the case uses:

```text
linear Kelvin-Voigt airway resistance
nonlinear Kelvin-Voigt airway resistance
airway inertia
```

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

Step 9.1 added direct-vs-generic structured assembly coverage for nonlinear Kelvin-Voigt airways. Existing Kelvin-Voigt airway sparse/tree comparison tests also continue to exercise the linear Kelvin-Voigt path.
