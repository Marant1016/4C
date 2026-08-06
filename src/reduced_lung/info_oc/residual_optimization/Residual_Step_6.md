# Residual Optimization Step 6

Date: 2026-08-06

## Scope

Implemented Step 6 from `residual_evaluation_optimization_guide.md`: reprofile after scalar residual optimizations before adding SIMD. No source code was changed for this step.

## Implemented Changes

- Rebuilt the release `4C` target so the scalar residual optimizations from earlier steps were included in the profile run.
- Reran the gen16 NewtonTree profile with `FOUR_C_REDUCED_LUNG_TREE_PROFILE=1`.
- Did not add SIMD kernels because Step 6 requires profiling before deciding whether SIMD work is justified.

## Profile Command

```bash
FOUR_C_REDUCED_LUNG_TREE_PROFILE=1 ./build/release/4C \
  /scratch/Rodriguez/workspace/4C/files/reduced_lung_lung_tree_gen16_500steps_newton_tree.4C.yaml \
  /tmp/opencode/residual_step6_profile/gen16_newton_tree_profile
```

## Profile Result

Key results from the run:

```text
newton_total_s: 9.11321
state_sync_s: 0.839854
residual_s: 1.23741
residual_clear_s: 0.241188
residual_airways_s: 0.347514
residual_terminal_units_s: 0.254647
residual_junctions_s: 0.531678
residual_boundary_conditions_s: 0.101938
residual_other_s: 0
residual_norm_s: 0.424793
residual_evaluations: 1500
tree_assembly_s: 1.37777
linear_solve_s: 4.43183
tree_solve_s: 4.43081
Calculation: 18.92
```

## Interpretation

- `residual_s` dropped substantially compared with the pre-residual-optimization baseline of about `5.94 s`.
- The largest residual callback phase is now `residual_junctions_s`.
- The next residual callback candidates by time are `residual_airways_s`, then `residual_terminal_units_s`, then `residual_boundary_conditions_s`.
- `residual_norm_s` and `residual_clear_s` are outside `residual_s` and are now visible enough to track, but Step 6 does not authorize changing them.
- The total Newton time is now dominated by `linear_solve_s` / `tree_solve_s`, not residual evaluation.

## SIMD Decision Point

Based on this profile, SIMD work should not start blindly. If SIMD is pursued next, the candidates should be considered in this order:

1. Bifurcation/connection junction residuals.
2. Rigid linear airway residuals.
3. Linear Kelvin-Voigt terminal-unit residuals.
4. Pressure boundary residuals.

The guide requirement remains: use SoA model arrays, padded or masked tails, and scalar fallback paths.

## Diff Summary

Focused Step 6 source diff:

```text
No source code changes.
```

Step 6 added this documentation file only.

## Verification

Build:

```bash
cmake --build build/release --target 4C --parallel 4
```

Result: passed.

Runtime profile:

```bash
FOUR_C_REDUCED_LUNG_TREE_PROFILE=1 ./build/release/4C \
  /scratch/Rodriguez/workspace/4C/files/reduced_lung_lung_tree_gen16_500steps_newton_tree.4C.yaml \
  /tmp/opencode/residual_step6_profile/gen16_newton_tree_profile
```

Result: passed, processor finished normally.

Tests: not run for this step because Step 6 made no source-code changes.

## Notes

- The profile output was written under `/tmp/opencode` to avoid inspecting workspace output directories.
- No build, output, log, or generated files were read to create this note.
