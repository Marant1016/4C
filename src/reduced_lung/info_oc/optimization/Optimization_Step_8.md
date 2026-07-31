# Optimization Step 8 - Vectorize Serial Top-Down Batches

## Scope

This step implements Phase 8 from the serial tree solver vectorization guide.

Only the serial `TreeNewtonLinearSolver` top-down recovery path was changed. The distributed `DistributedTreeNewtonLinearSolver` implementation was intentionally left unchanged.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Added Top-Down Batch Workspace

The serial solver now owns reusable lane buffers for grouped top-down recovery:

```cpp
std::vector<double> top_down_inlet_pressure_;
std::vector<double> top_down_unknown_values_;
std::vector<double> top_down_outlet_pressure_;
std::vector<double> top_down_child_pressure_;
```

These buffers are sized during symbolic setup from the largest top-down group.

### Replaced Scalar Top-Down Element Loop

Top-down recovery now runs through `recover_top_down_group(...)` for each top-down layer group.

The grouped recovery performs these phases:

```text
gather inlet-pressure corrections for all group lanes
write inlet-pressure delta entries
for each unknown slot:
  compute correction values for all group lanes
  scatter correction values to delta
for non-leaf groups:
  compute outlet pressure for all group lanes
  for each child slot:
    compute child inlet pressure for all group lanes
    stamp/write child inlet pressures
```

The correction formula is unchanged:

```cpp
value = slope * inlet_pressure + intercept;
```

### Preserved Stamp-Based Safety

The inlet-pressure stamp checks from Step 4 are still used. A group lane cannot read an inlet-pressure correction unless it was written during the current solve stamp.

Child inlet-pressure writes still use the same duplicate-write assertion.

## Behavior Preserved

- The top-down mathematical recovery formula is unchanged.
- Top-down layer dependency order is unchanged.
- Elements are only processed in the existing same-layer groups.
- Delta scatter still uses `replace_global_value(...)`.
- Bottom-up batch behavior from Step 7 is unchanged.
- Coefficient access behavior from Step 6 is unchanged.
- No SIMD intrinsics were added.
- The distributed tree solver was not changed.

## Guide Status

The serial vectorization guide phases are now implemented in the serial solver path:

- Phase 1: baseline measurements documented in `Optimization_Step_1.md`.
- Phase 2: serial plan/workspace SoA layout.
- Phase 3: bottom-up and top-down layer groups.
- Phase 4: avoidable serial hot-path clears removed.
- Phase 5: `2x2` and `3x3` batch dense solvers with scalar fallback.
- Phase 6: virtual coefficient-provider dispatch removed from serial hot loops.
- Phase 7: bottom-up traversal uses group assemble, batch solve, and group relation write phases.
- Phase 8: top-down traversal uses grouped lane-wise recovery phases.

The remaining structured-coefficient row-entry search in `TreeLinearization::value(...)` is still a follow-up optimization. The guide listed precomputed coefficient locations and direct tree block storage as better follow-ups rather than the minimum serial vectorization milestone.

## Intended Benefit

This step completes the grouped serial tree traversal design: both bottom-up condensation and top-down recovery now operate by layer groups instead of fully scalar per-element traversal. The loops are still explicit C++ loops, but their structure is regular and same-shape, making them more suitable for compiler vectorization and later targeted micro-optimizations.

## Verification

Verification passed:

```text
git diff --check
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

## Post-Optimization Benchmark Rerun

The release benchmark target was rebuilt and the Step 1 baseline benchmark commands were rerun on 2026-07-31 after Steps 2-8:

```text
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/release --target benchmarktests_reduced_lung --parallel 4
./build/release/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/LinearSolve/BalancedAirways/StructuredTree --benchmark_repetitions=5 --benchmark_min_time=0.01s
./build/release/tests/benchmarktests_reduced_lung --benchmark_filter="ReducedLung/FullSolve/(SingleTerminalUnit|SerialAirways|BalancedAirways)/NewtonTree" --benchmark_repetitions=5 --benchmark_min_time=0.01s
```

The benchmark run emitted the same expected authorization, CPU-scaling, and PHG redistribution warnings as the baseline. The structured linear-solve run also had high CV for several cases (`StructuredTree/3` through `StructuredTree/5`), so these values should be treated as directional.

The full-solve raw output was saved by the tool harness at:

```text
/data/home/Rodriguez/.local/share/opencode/tool-output/tool_fb7aa8934001yl8lmNF58Z32zi
```

### Structured Tree Linear-Solve Means After Steps 2-8

| Benchmark | Elements | Dofs | Time | tree_solve_s | tree_bottom_up_s | tree_top_down_s | tree_dense_s | tree_lookup_s | tree_dense_solves | tree_lookups | tree_workspace_dofs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `StructuredTree/2_mean` | 3 | 9 | 3.62 us | 3.57132 us | 2.57131 us | 673.297 ns | 144.538 ns | 1.13985 us | 3 | 25 | 6 |
| `StructuredTree/3_mean` | 7 | 21 | 7.32 us | 7.28161 us | 5.68451 us | 1.30619 us | 238.468 ns | 2.44153 us | 7 | 61 | 14 |
| `StructuredTree/4_mean` | 15 | 45 | 12.2 us | 12.2106 us | 9.91442 us | 2.06206 us | 346.583 ns | 4.19262 us | 15 | 133 | 30 |
| `StructuredTree/5_mean` | 31 | 93 | 20.6 us | 20.5871 us | 17.1735 us | 3.21247 us | 538.696 ns | 7.17352 us | 31 | 277 | 62 |

### Structured Tree Linear-Solve Delta Versus Step 1

| Benchmark | Time | tree_solve_s | tree_bottom_up_s | tree_top_down_s | tree_dense_s | tree_lookup_s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `StructuredTree/2_mean` | +91.5% | +91.0% | +76.7% | +174.6% | +6.3% | +114.9% |
| `StructuredTree/3_mean` | +23.9% | +24.0% | +19.1% | +51.9% | -36.0% | +37.5% |
| `StructuredTree/4_mean` | +9.9% | +10.1% | +5.9% | +36.8% | -51.6% | +22.8% |
| `StructuredTree/5_mean` | -0.5% | -0.2% | -3.6% | +24.7% | -64.0% | +10.4% |

### Full `NewtonTree` Solve Means After Steps 2-8

| Benchmark | Elements | Dofs | Time | newton_total_s | linear_solve_s | tree_solve_s | tree_bottom_up_s | tree_top_down_s | tree_dense_s | tree_lookup_s | tree_assembly_s | residual_s | state_sync_s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `SingleTerminalUnit/NewtonTree_mean` | 1 | 3 | 17.4 us | 14.2085 us | 2.77243 us | 2.57056 us | 1.48741 us | 452.754 ns | 352.118 ns | 557.193 ns | 1.0651 us | 2.68723 us | 6.53161 us |
| `SerialAirways/NewtonTree_mean` | 9 | 27 | 43.5 us | 40.3267 us | 14.9048 us | 14.7076 us | 11.739 us | 2.32664 us | 1.21582 us | 4.26844 us | 4.45596 us | 5.04992 us | 14.6104 us |
| `BalancedAirways/NewtonTree_mean` | 15 | 45 | 57.3 us | 54.3351 us | 21.7832 us | 21.6116 us | 17.307 us | 3.67703 us | 904.945 ns | 7.02548 us | 6.16467 us | 5.89611 us | 19.1609 us |

All three full-solve benchmark cases still reported `nonlinear_iterations=2`.

### Full `NewtonTree` Solve Delta Versus Step 1

| Benchmark | Time | newton_total_s | linear_solve_s | tree_solve_s | tree_bottom_up_s | tree_top_down_s | tree_dense_s | tree_lookup_s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `SingleTerminalUnit/NewtonTree_mean` | +5.5% | +4.9% | +8.6% | +5.6% | -1.7% | +50.1% | -2.0% | -9.3% |
| `SerialAirways/NewtonTree_mean` | +2.6% | +2.1% | +0.9% | +0.5% | -2.7% | +23.0% | -3.5% | -1.6% |
| `BalancedAirways/NewtonTree_mean` | -22.4% | -23.3% | -25.3% | -25.5% | -28.1% | -11.1% | -53.2% | -22.2% |

### Post-Optimization Conclusion

The full Step 2-8 series did not produce uniform speedups in this noisy benchmark environment. The largest structured linear-solve case (`StructuredTree/5_mean`) is essentially flat against Step 1 (`20.5871 us` tree solve after Steps 2-8 versus `20.6226 us` baseline), with lower bottom-up and dense-solve time offset by higher top-down and lookup time.

The dense batch solver work is visible on larger linear cases: `tree_dense_s` dropped by about `52%` for `StructuredTree/4_mean` and `64%` for `StructuredTree/5_mean`. Smaller structured-tree cases regressed because the added grouping/staging overhead dominates their tiny work sizes.

The full `BalancedAirways/NewtonTree_mean` case improved substantially (`tree_solve_s` down about `25.5%`, `newton_total_s` down about `23.3%`), but the current run's CPU-scaling warning and `12.17%` CV mean this should be treated as directional rather than final proof. The single-terminal and serial-airways full-solve cases are effectively flat to slightly slower.

The next useful optimization should target the remaining structured-coefficient lookup overhead and top-down staging overhead. The most direct follow-up is still the guide's precomputed coefficient-location/direct tree-block storage work, because `tree_lookup_s` remains a large share of the measured tree solve and increased in the standalone structured linear-solve rerun.
