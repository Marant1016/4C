# Residual Optimization Step 1

Date: 2026-08-06

## Scope

Implemented Step 1 from `residual_evaluation_optimization_guide.md`: residual phase profiling only. No residual algorithm, write path, allocation strategy, or numerical behavior was changed.

## Implemented Changes

- Added residual profiling fields to `NewtonSolverProfile` for clear time, per-phase residual callback time, norm time, and residual evaluation count.
- Added `NamedResidualAssembler` to `ReducedLungAssemblyPipeline`.
- Registered named residual phases in the default reduced-lung assembly pipeline:
  - airways;
  - terminal units;
  - junctions;
  - boundary conditions.
- Updated `NewtonSolver::assemble_residual_for_current_state()` to time:
  - `residual_.put_scalar(0.0)` as `residual_clear_s`;
  - named residual callback groups as phase-specific residual timers;
  - unlabelled fallback residual callbacks as `residual_other_s`;
  - residual callback total as the existing `residual_s`;
  - `residual_.norm_2(...)` as `residual_norm_s`.
- Guarded the new fine-grained timer starts behind the profile pointer so non-profiled solves do not pay extra per-phase timing overhead.
- Updated the reduced-lung tree profile summary to print the new residual timing fields.

## Modified Files

- `src/reduced_lung/src/4C_reduced_lung_solver_profile.hpp`
- `src/reduced_lung/src/4C_reduced_lung_helpers.hpp`
- `src/reduced_lung/src/4C_reduced_lung_helpers.cpp`
- `src/reduced_lung/src/4C_reduced_lung_newton_solver.cpp`
- `src/reduced_lung/src/4C_reduced_lung_main.cpp`

## Diff Summary

From the focused git diff:

```text
5 files changed, 96 insertions(+), 3 deletions(-)
```

Change distribution:

- `4C_reduced_lung_solver_profile.hpp`: new Newton residual timing counters.
- `4C_reduced_lung_helpers.hpp`: named residual assembler type and vector.
- `4C_reduced_lung_helpers.cpp`: phase labels attached to default residual assemblers.
- `4C_reduced_lung_newton_solver.cpp`: residual clear, phase, total callback, fallback, norm, and evaluation-count timing.
- `4C_reduced_lung_main.cpp`: tree-profile printout extended with residual phase counters.

## Verification

Build:

```bash
cmake --build build/debug --target unittests_reduced_lung --parallel 4
```

Result: passed.

Tests:

```bash
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

Result: both passed.

## Notes

- The existing `residual_s` remains the total time spent in residual assembler callbacks.
- `residual_clear_s` and `residual_norm_s` are separate from `residual_s`.
- Existing unlabelled residual callbacks are still supported and are counted under `residual_other_s` when no named residual assemblers are present.
- No residual evaluation algorithm was optimized in this step; this prepares the measurements needed for later steps.
