# Implementation NextStep 7 - Default Solver Policy

## Scope

This step formalizes the reduced-lung nonlinear solver default policy without changing any YAML input
files.

Policy retained:

- `Nox` remains the default nonlinear solver workflow.
- `NewtonSparse` remains an explicit opt-in reference/debugging workflow.
- `NewtonTree` remains an explicit opt-in tree-solver workflow.
- Existing YAML files that omit `nonlinear_solver` continue to run with `Nox`.
- YAML variants that set `nonlinear_solver: NewtonSparse` or `nonlinear_solver: NewtonTree` continue
  to select those workflows explicitly.

## Files Modified

- `src/reduced_lung/src/4C_reduced_lung_input.hpp`
- `src/reduced_lung/tests/4C_reduced_lung_input_pipeline_test.cpp`

## Files Added

- `src/reduced_lung/info_oc/Implementation_NextStep_7.md`

## YAML Files

No YAML input files were changed.

The intended input behavior is:

```yaml
REDUCED DIMENSIONAL LUNG DYNAMIC:
  # nonlinear_solver omitted -> Nox
```

Explicit tree opt-in remains:

```yaml
REDUCED DIMENSIONAL LUNG DYNAMIC:
  nonlinear_solver: NewtonTree
```

Explicit sparse Newton opt-in remains:

```yaml
REDUCED DIMENSIONAL LUNG DYNAMIC:
  nonlinear_solver: NewtonSparse
```

## Code Changes

`ReducedLungParameters::Dynamics::nonlinear_solver` still defaults to:

```cpp
NonlinearSolverType::Nox
```

A short comment was added beside the default to document that this is intentional and that
`NewtonSparse`/`NewtonTree` are explicit opt-in workflows.

No runtime solver-selection logic was changed.

## Test Changes

Extended reduced-lung input-policy coverage in:

```text
src/reduced_lung/tests/4C_reduced_lung_input_pipeline_test.cpp
```

Existing coverage retained:

```text
ReducedLungInputPipelineTest.NonlinearSolverDefaultsToNox
```

New coverage added:

```text
ReducedLungInputPipelineTest.NewtonNonlinearSolversAreExplicitOptInValues
```

The new test verifies that `NewtonSparse` and `NewtonTree` remain valid explicit opt-in enum values
without changing the default.

## Rationale

`NewtonTree` is now validated across more serial and distributed correction paths, but `Nox` should
remain the safe default until the tree path has broader representative input coverage and performance
evidence.

Keeping `Nox` as default preserves baseline behavior for existing input decks. Keeping the Newton
paths explicit makes solver comparisons and debugging intentional and reproducible.

## Current Policy Summary

Recommended usage:

- Omit `nonlinear_solver` for baseline `Nox` behavior.
- Use `nonlinear_solver: NewtonSparse` when a custom Newton sparse reference is needed.
- Use `nonlinear_solver: NewtonTree` when testing or benchmarking the structured tree solver.
- Do not remove NOX unless the replacement path fully covers required functionality and an explicit
  fallback strategy exists.
