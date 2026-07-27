# Optimization Step 1 - Serial Tree Solver Layout Cleanup

## Scope

This step applies the first low-risk optimization pass to the serial reduced-lung tree linear solver.

Only the serial `TreeNewtonLinearSolver` path was changed. The distributed `DistributedTreeNewtonLinearSolver` path was intentionally left unchanged because it is a separate implementation with its own communication and workspace logic.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Flattened Local Dense Matrix Storage

The serial `ElementWorkspace::matrix` storage changed from nested vectors:

```cpp
std::vector<std::vector<double>> matrix;
```

to one contiguous flat vector:

```cpp
std::vector<double> matrix;
```

The dense solver now indexes entries as:

```cpp
matrix[row * block_size + col]
```

This removes one layer of pointer indirection per matrix row and improves locality for the small dense local systems used during bottom-up condensation.

### Fixed Serial Child-Interface Storage

The serial tree solver supports at most two children per parent. The per-element child-interface storage changed from a dynamically allocated vector to fixed-size storage:

```cpp
std::array<ChildInterfacePlan, 2> child_interfaces;
int child_interface_count;
```

The matching numeric child-pressure recovery arrays in `ElementWorkspace` were also changed to fixed-size arrays:

```cpp
std::array<double, 2> child_pressure_slope;
std::array<double, 2> child_pressure_intercept;
```

This avoids per-element heap storage for child-interface metadata in the serial path.

### Split Serial Subtree Relations

The serial subtree relation storage changed from an Array-of-Structs layout:

```cpp
std::vector<SubtreeRelation> subtree_relations;
```

to separate arrays:

```cpp
std::vector<double> subtree_relation_G;
std::vector<double> subtree_relation_h;
```

This is a small Struct-of-Arrays step and prepares the bottom-up loop for later layer-wise batching or vectorization.

## Behavior Preserved

- The bottom-up/top-down tree algorithm was not changed.
- The Newton correction convention remains `J * delta = -F`.
- The structured and sparse coefficient-source behavior is unchanged.
- The serial tree solver still reuses its symbolic plan and per-element workspaces.
- The distributed tree solver implementation was not changed.

## Intended Benefit

This step reduces pointer chasing and dynamic per-element storage in the serial solver hot path. It is intended as preparation for later optimizations such as specialized `2x2`/`3x3` local solvers and direct structured coefficient access.

## Verification And Baseline Measurement

Verification run successfully:

```text
git diff --check
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
```

The serial reduced-lung unit-test target passed.

Additional Phase 1 baseline measurements were run on 2026-07-27 from repository state:

```text
932b3472c0783174ed1b5d66ee7f751e16b0d6a3
```

Before updating this note, `git status --short` produced no output.

The release reduced-lung benchmark target was built with:

```text
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/release --target benchmarktests_reduced_lung --parallel 4
```

Serial structured tree linear-solve baseline command:

```text
./build/release/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/LinearSolve/BalancedAirways/StructuredTree --benchmark_repetitions=5 --benchmark_min_time=0.01s
```

Serial `NewtonTree` full-solve baseline command:

```text
./build/release/tests/benchmarktests_reduced_lung --benchmark_filter="ReducedLung/FullSolve/(SingleTerminalUnit|SerialAirways|BalancedAirways)/NewtonTree" --benchmark_repetitions=5 --benchmark_min_time=0.01s
```

The benchmark run emitted expected CPU-scaling and PHG redistribution warnings, so the values should be treated as baseline guidance rather than final performance numbers.

### Structured Tree Linear-Solve Means

| Benchmark | Elements | Dofs | Time | tree_solve_s | tree_bottom_up_s | tree_top_down_s | tree_dense_s | tree_lookup_s | tree_dense_solves | tree_lookups | tree_workspace_dofs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `StructuredTree/2_mean` | 3 | 9 | 1.89 us | 1.86968 us | 1.45522 us | 245.232 ns | 135.933 ns | 530.609 ns | 3 | 25 | 6 |
| `StructuredTree/3_mean` | 7 | 21 | 5.91 us | 5.87335 us | 4.7734 us | 859.658 ns | 372.418 ns | 1.77573 us | 7 | 61 | 14 |
| `StructuredTree/4_mean` | 15 | 45 | 11.1 us | 11.0916 us | 9.36255 us | 1.50767 us | 716.339 ns | 3.4132 us | 15 | 133 | 30 |
| `StructuredTree/5_mean` | 31 | 93 | 20.7 us | 20.6226 us | 17.8217 us | 2.57538 us | 1.49776 us | 6.49799 us | 31 | 277 | 62 |

### Full `NewtonTree` Solve Means

| Benchmark | Elements | Dofs | Time | newton_total_s | linear_solve_s | tree_solve_s | tree_bottom_up_s | tree_top_down_s | tree_dense_s | tree_lookup_s | tree_assembly_s | residual_s | state_sync_s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `SingleTerminalUnit/NewtonTree_mean` | 1 | 3 | 16.5 us | 13.5473 us | 2.55404 us | 2.43342 us | 1.51263 us | 301.644 ns | 359.465 ns | 614.491 ns | 1.09612 us | 2.28722 us | 6.37775 us |
| `SerialAirways/NewtonTree_mean` | 9 | 27 | 42.4 us | 39.5042 us | 14.7652 us | 14.6393 us | 12.0637 us | 1.89213 us | 1.26031 us | 4.3381 us | 4.57375 us | 4.5209 us | 14.2906 us |
| `BalancedAirways/NewtonTree_mean` | 15 | 45 | 73.8 us | 70.8315 us | 29.1507 us | 29.0119 us | 24.0599 us | 4.13465 us | 1.93199 us | 9.02792 us | 7.95481 us | 7.3252 us | 24.8369 us |

All three full-solve benchmark cases reported `nonlinear_iterations=2`.

### Baseline Conclusion

For the serial structured tree linear solve, bottom-up condensation is the dominant phase. In the largest measured linear-solve case (`StructuredTree/5_mean`), `tree_bottom_up_s` was about `17.8217 us` of `20.6226 us` total tree solve time. Coefficient lookup was also significant at `6.49799 us`, while dense local solves were only `1.49776 us`.

For the full `BalancedAirways/NewtonTree_mean` solve, the tree solve took `29.0119 us` of `70.8315 us` custom Newton time. State synchronization (`24.8369 us`), structured tree assembly (`7.95481 us`), and residual assembly (`7.3252 us`) are also visible costs.

The next optimization phase should therefore prioritize SoA/grouped bottom-up traversal and reduced coefficient lookup overhead before focusing only on dense local solver kernels.
