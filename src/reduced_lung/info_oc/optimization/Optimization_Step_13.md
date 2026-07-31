# Optimization Step 13 - Select Scalar Versus Grouped Solver By Tree Size

## Scope

This step adds an internal serial tree solver policy that chooses scalar full-tree traversal for small trees and keeps the grouped/vectorization-friendly traversal for larger trees.

Only the serial `TreeNewtonLinearSolver` path was changed. The distributed `DistributedTreeNewtonLinearSolver` implementation was intentionally left unchanged.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Added Tree-Size Policy

The serial solver now stores an internal traversal policy flag:

```cpp
bool use_scalar_tree_solve_ = false;
```

`build_symbolic_plan()` sets this flag from the number of tree elements:

```cpp
constexpr int scalar_tree_element_threshold = 7;
use_scalar_tree_solve_ = element_count <= scalar_tree_element_threshold;
```

Trees with up to seven elements use scalar bottom-up and top-down traversal. Larger trees keep the grouped traversal from the previous steps.

### Added Scalar Full-Tree Bottom-Up Traversal

Small trees now traverse `tree_metadata_.bottom_up_layers` directly. For each element, the solver:

```text
assembles the element rows
applies child pressure/flow condensation if needed
solves the local dense block with the scalar dense solver
writes the subtree relation
```

The grouped path is still used for larger trees, including the Step 12 specialized `2x2` assembly helpers and the `2x2`/`3x3` batch dense solvers.

### Reused Scalar Top-Down Recovery For Small Trees

Small trees now traverse `tree_metadata_.top_down_layers` directly and call the Step 11 scalar `recover_top_down_element(...)` path for each element.

Larger trees keep the grouped top-down path, including the Step 11 small-group scalar fallback inside grouped traversal.

## Behavior Preserved

- The bottom-up condensation mathematics is unchanged.
- The top-down recovery formula is unchanged.
- Coefficient storage and direct coefficient values from Steps 9-10 are unchanged.
- Required coefficient assertions are unchanged.
- Sparse-Jacobian fallback behavior is unchanged.
- The distributed tree solver was not changed.

## Verification

Verification passed:

```text
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

The reduced-lung unit test and focused `newton_tree` input tests passed.

## Release Benchmark Rerun

The release benchmark target was rebuilt and the same benchmark commands from `Optimization_Step_8.md` were rerun on 2026-07-31:

```text
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/release --target benchmarktests_reduced_lung --parallel 4
./build/release/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/LinearSolve/BalancedAirways/StructuredTree --benchmark_repetitions=5 --benchmark_min_time=0.01s
./build/release/tests/benchmarktests_reduced_lung --benchmark_filter="ReducedLung/FullSolve/(SingleTerminalUnit|SerialAirways|BalancedAirways)/NewtonTree" --benchmark_repetitions=5 --benchmark_min_time=0.01s
```

The benchmark run emitted the same expected authorization, CPU-scaling, and PHG redistribution warnings as the Step 8 rerun. Several structured linear-solve cases still had high CV, so the values should be treated as directional.

The full-solve raw output was saved by the tool harness at:

```text
/data/home/Rodriguez/.local/share/opencode/tool-output/tool_fb82b199b001y530wJu3P6Cfev
```

### Structured Tree Linear-Solve Means After Step 13

| Benchmark | Elements | Dofs | Time | tree_solve_s | tree_bottom_up_s | tree_top_down_s | tree_dense_s | tree_lookup_s | tree_dense_solves | tree_lookups | tree_workspace_dofs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `StructuredTree/2_mean` | 3 | 9 | 0.905 us | 0.875661 us | 0.426304 us | 0.3121 us | 0.144884 us | 0 | 3 | 0 | 6 |
| `StructuredTree/3_mean` | 7 | 21 | 1.61 us | 1.58314 us | 0.920437 us | 0.550563 us | 0.315173 us | 0 | 7 | 0 | 14 |
| `StructuredTree/4_mean` | 15 | 45 | 2.64 us | 2.61331 us | 1.31004 us | 1.19189 us | 0.290184 us | 0 | 15 | 0 | 30 |
| `StructuredTree/5_mean` | 31 | 93 | 5.96 us | 5.93419 us | 2.92978 us | 2.86843 us | 0.524934 us | 0 | 31 | 0 | 62 |

### Structured Tree Linear-Solve Delta Versus Step 8

| Benchmark | Time | tree_solve_s | tree_bottom_up_s | tree_top_down_s | tree_dense_s | tree_lookup_s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `StructuredTree/2_mean` | -75.0% | -75.5% | -83.4% | -53.6% | +0.2% | -100.0% |
| `StructuredTree/3_mean` | -78.0% | -78.3% | -83.8% | -57.9% | +32.2% | -100.0% |
| `StructuredTree/4_mean` | -78.4% | -78.6% | -86.8% | -42.2% | -16.3% | -100.0% |
| `StructuredTree/5_mean` | -71.1% | -71.2% | -82.9% | -10.7% | -2.6% | -100.0% |

### Full `NewtonTree` Solve Means After Step 13

| Benchmark | Elements | Dofs | Time | newton_total_s | linear_solve_s | tree_solve_s | tree_bottom_up_s | tree_top_down_s | tree_dense_s | tree_lookup_s | tree_assembly_s | residual_s | state_sync_s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `SingleTerminalUnit/NewtonTree_mean` | 1 | 3 | 17.6 us | 14.3675 us | 1.80771 us | 1.61999 us | 0.722591 us | 0.372173 us | 0.373452 us | 0 | 1.65072 us | 2.67415 us | 7.01508 us |
| `SerialAirways/NewtonTree_mean` | 9 | 27 | 37.5 us | 34.1278 us | 6.55466 us | 6.36901 us | 3.74432 us | 2.11267 us | 1.29635 us | 0 | 6.19584 us | 4.79051 us | 15.2343 us |
| `BalancedAirways/NewtonTree_mean` | 15 | 45 | 51.6 us | 48.4471 us | 9.06519 us | 8.89066 us | 4.31666 us | 4.00306 us | 0.976806 us | 0 | 9.62805 us | 6.4071 us | 21.9213 us |

All three full-solve benchmark cases still reported `nonlinear_iterations=2`.

### Full `NewtonTree` Solve Delta Versus Step 8

| Benchmark | Time | newton_total_s | linear_solve_s | tree_solve_s | tree_bottom_up_s | tree_top_down_s | tree_dense_s | tree_lookup_s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `SingleTerminalUnit/NewtonTree_mean` | +1.1% | +1.1% | -34.8% | -37.0% | -51.4% | -17.8% | +6.1% | -100.0% |
| `SerialAirways/NewtonTree_mean` | -13.8% | -15.4% | -56.0% | -56.7% | -68.1% | -9.2% | +6.6% | -100.0% |
| `BalancedAirways/NewtonTree_mean` | -9.9% | -10.8% | -58.4% | -58.9% | -75.1% | +8.9% | +7.9% | -100.0% |

### Benchmark Conclusion

The Step 9-13 follow-up series materially improved solve-time tree counters versus the Step 8 rerun. The largest standalone structured linear case (`StructuredTree/5_mean`) dropped from `20.5871 us` tree solve time in Step 8 to `5.93419 us` after Step 13. The smaller scalar-policy cases also improved substantially: `StructuredTree/2_mean` dropped from `3.57132 us` to `0.875661 us`, and `StructuredTree/3_mean` dropped from `7.28161 us` to `1.58314 us`.

Full-solve tree time also improved: `BalancedAirways/NewtonTree_mean` tree solve time dropped from `21.6116 us` in Step 8 to `8.89066 us` after Step 13. The total custom Newton time improved less because some cost moved into structured tree assembly and other full-solve phases. This is visible in `tree_assembly_s`, which rose for the full solves after direct coefficient storage.

`tree_lookup_s` and `tree_lookups` are now zero for the structured serial benchmarks because the solve loop reads direct coefficient values instead of calling row-entry lookup helpers. That is expected for Steps 10-13, but it means the cost moved earlier into structured coefficient preparation rather than disappearing entirely.

Compared with the Step 1 baseline, the final `BalancedAirways/NewtonTree_mean` tree solve time is down from `29.0119 us` to `8.89066 us`, and full custom Newton time is down from `70.8315 us` to `48.4471 us`. The single-terminal full solve remains slightly slower in total wall time than Step 8 despite lower tree solve time, so the scalar threshold should be revisited if single-terminal total runtime becomes the primary target.
