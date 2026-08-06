# Assemply Step 9.1

## Goal

Add correctness coverage before optimizing the remaining nonlinear and non-gen16 model coefficient loops.

This step changes tests only. It does not change model formulas or production assembly behavior.

## Why

The Step 9 plan targets model paths that were not the main gen16 case:

```text
nonlinear rigid airways
nonlinear Kelvin-Voigt airways
Ogden terminal units
Four-element Maxwell terminal units
Four-element Maxwell + Ogden terminal units
```

Before optimizing those loops, the test suite should prove that the direct structured tree assembly path gives the same coefficients as the generic `TreeLinearization` path.

## Test Changes

File:

```text
src/reduced_lung/tests/4C_reduced_lung_tree_linear_solver_test.cpp
```

Added fixture helpers:

```text
make_nonlinear_kelvin_voigt_airway_parameters(...)
make_ogden_terminal_unit_parameters(...)
make_four_element_maxwell_ogden_terminal_unit_parameters(...)
```

Reused existing helpers:

```text
make_nonlinear_airway_parameters(...)
make_four_element_maxwell_terminal_unit_parameters(...)
compare_direct_and_generic_structured_coefficients(...)
```

## New Tests

Added direct-vs-generic structured assembly tests:

```text
DirectStructuredAssemblyNonlinearRigidAirwaysMatchesGenericPath
DirectStructuredAssemblyNonlinearKelvinVoigtAirwaysMatchesGenericPath
DirectStructuredAssemblyOgdenTerminalUnitMatchesGenericPath
DirectStructuredAssemblyFourElementMaxwellTerminalUnitMatchesGenericPath
DirectStructuredAssemblyFourElementMaxwellOgdenTerminalUnitMatchesGenericPath
```

Each test checks:

```text
1. Generic TreeLinearization assembly.
2. Direct TreeNewtonLinearSolver coefficient assembly.
3. Every coefficient consumed by TreeNewtonLinearSolver matches.
4. Direct structured tree correction matches the sparse solver correction.
5. Direct structured path performs zero coefficient lookups.
```

The nonlinear tests use seeded nonzero states so the nonlinear formulas are actually exercised.

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

## Result

The direct structured assembly path is now covered for the remaining model families that Step 9 intends to optimize. Future Step 9 production changes can be checked against these tests before relying on performance measurements.
