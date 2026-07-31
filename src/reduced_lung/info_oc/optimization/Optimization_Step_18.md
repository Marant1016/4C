# Optimization Step 18 - Gate Profiling Timers And Keep Thresholds Stable

## Scope

This step removes unnecessary timer calls from the serial `TreeNewtonLinearSolver` path when no tree profile is attached.

Only the serial tree solver implementation was changed. The distributed tree solver was intentionally left unchanged.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Gated Serial Tree Solver Timer Starts

The serial solver now starts detailed timers only when `profile_ != nullptr`.

Before this step, these timers called `Clock::now()` unconditionally:

- total tree solve timer
- scalar dense-solve timer
- bottom-up timer
- grouped `2x2` dense-solve timer
- grouped `3x3` dense-solve timer
- top-down timer

They now use the guarded form:

```cpp
const auto timer_start = profile_ != nullptr ? Clock::now() : Clock::time_point{};
```

The elapsed-time update is still inside the existing `if (profile_ != nullptr)` blocks, so profile counters keep the same behavior when profiling is enabled.

### Thresholds Left Unchanged

No scalar/grouped threshold changes were made in this step.

Current values remain:

```cpp
constexpr int scalar_tree_element_threshold = 7;
constexpr int top_down_scalar_group_threshold = 2;
```

Reason: retuning these thresholds should be based on a release benchmark sweep after Steps 14-17. No release benchmark sweep was run for Step 18, so changing thresholds would be guesswork.

## Behavior Preserved

- The Newton equations are unchanged.
- The bottom-up and top-down solver mathematics are unchanged.
- Detailed benchmark/profile counters are preserved when `profile_` is attached.
- Production/no-profile serial solves avoid the gated timer calls.
- Scalar/grouped solver thresholds are unchanged.
- The distributed tree solver was not changed.

## Expected Performance Effect

This step targets production/no-profile serial `NewtonTree` runs by avoiding unnecessary `Clock::now()` calls.

Benchmarks that attach a `TreeNewtonLinearSolverProfile` should keep the same detailed timing counters and still pay the timing cost intentionally.

## Verification

Verification passed:

```text
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
git diff --check
```

The `ctest` commands were run from `build/debug`.

Results:

- `unittests_reduced_lung`: passed, `1/1` tests.
- Focused `newton_tree` input tests: passed, `7/7` tests.
- `git diff --check`: passed.

## Benchmark Status

Release benchmarks were not run for this step.

When benchmarking later, compare against Step 17 and Step 13 with focus on:

- no-profile production runtime if practical
- `tree_solve_s`
- `tree_bottom_up_s`
- `tree_top_down_s`
- `tree_dense_s`
- `newton_total_s`
- nonlinear iteration count

If a future threshold sweep is run, test at least:

- `scalar_tree_element_threshold`: `3`, `7`, `15`
- `top_down_scalar_group_threshold`: `1`, `2`, `4`

Keep the current thresholds unless a release benchmark sweep shows a clear total-runtime improvement.

## Conclusion

Step 18 keeps solver behavior and thresholds unchanged while removing detailed timer overhead from no-profile serial tree solves. Threshold retuning remains a benchmark-driven follow-up rather than a guessed code change.
