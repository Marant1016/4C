# Implementation NextStep 8 - Distributed Tree Benchmark Input And Harness

## Scope

This step adds only new files. No production solver code was modified.

The goal is to make it possible to exercise and benchmark the MPI-capable `NewtonTree` path that uses `DistributedTreeNewtonLinearSolver`.

## Files Added

### `src/reduced_lung/benchmark_tests/4C_reduced_lung_distributed_tree_solver_benchmark.cpp`

Added a new Google Benchmark source file for distributed reduced-lung tree-solver measurements.

Because `src/reduced_lung/benchmark_tests/CMakeLists.txt` uses `four_c_auto_define_benchmark_tests()`, the new `.cpp` file is picked up by the existing `benchmarktests_reduced_lung` target without modifying CMake.

The benchmark file builds distributed reduced-lung fixtures directly in C++. It does not read YAML input files.

Benchmark groups added:

```text
ReducedLung/Distributed/FullSolve/BalancedAirways/Nox
ReducedLung/Distributed/FullSolve/BalancedAirways/NewtonSparse
ReducedLung/Distributed/FullSolve/BalancedAirways/NewtonTree
ReducedLung/Distributed/LinearSolve/BalancedAirways/StructuredTree
```

The full-solve benchmarks use a 4-level balanced rigid-airway tree. The linear-solve benchmark uses balanced rigid-airway trees with 2, 3, 4, and 5 levels.

The physical model is intentionally simple:

- airway elements only,
- linear resistance,
- rigid walls,
- no inertia,
- pressure boundary condition at the root inlet,
- pressure boundary condition at every leaf outlet.

This matches the current serial benchmark style and focuses the measurement on solver and communication cost rather than model complexity.

## MPI Timing Design

The new distributed benchmarks use fixed iterations and Google Benchmark manual timing:

```text
UseManualTime()
Iterations(...)
```

Each measured iteration does:

```text
MPI_Barrier
start = MPI_Wtime()
run solve operation
MPI_Barrier
elapsed = MPI_Wtime() - start
MPI_Allreduce(elapsed, max_elapsed, MPI_MAX)
state.SetIterationTime(max_elapsed)
```

The `MPI_MAX` reduction reports distributed wall-clock time, i.e. the time of the slowest rank.

This avoids the main MPI benchmark problem where Google Benchmark could otherwise choose a different number of iterations on different ranks and deadlock inside collective communication.

## Recorded Counters

The benchmark records common counters:

```text
mpi_ranks
elements
dofs
equations
cross_rank_edges
nonlinear_iterations
final_residual
```

For NOX, it records profile counters such as residual evaluations, Jacobian evaluations, residual assembly time, sparse Jacobian assembly time, sparse matrix completion time, state synchronization time, and total solve time.

For `NewtonSparse`, it records custom Newton counters and sparse linear solve time.

For `NewtonTree`, it records custom Newton counters plus distributed tree-solver counters:

```text
tree_solve_s
tree_bottom_up_s
tree_top_down_s
tree_dense_s
tree_lookup_s
tree_comm_s
tree_dense_solves
tree_lookups
relation_messages
pressure_messages
scatter_messages
communication_bytes
tree_workspace_dofs
tree_max_block
```

Tree message counters correspond to the currently implemented distributed algorithm:

- bottom-up `(G, h)` boundary relation messages,
- top-down inlet-pressure messages,
- final correction scatter messages.

## Distributed Runtime YAML Input

### `tests/input_files/reduced_lung_distributed_bifurcation_newton_tree.4C.yaml`

Added an NP2 runtime input file for testing the distributed `NewtonTree` executable path.

The file uses the existing bifurcation field data:

```text
reduced_lung_aw_bifurcation_fields.json
```

and sets:

```yaml
nonlinear_solver: NewtonTree
```

It is registered in `tests/list_of_tests.cmake` as an NP2 CTest input.

Manual example run:

```text
mpirun -np 2 ./4C ../../tests/input_files/reduced_lung_distributed_bifurcation_newton_tree.4C.yaml ../../../output_np2/
```

With `comm_size > 1`, runtime solver construction should select:

```text
DistributedTreeNewtonLinearSolver
```

## How To Run The Distributed Benchmark

Build the existing benchmark target after configuring Google Benchmark support:

```text
cmake --build build/debug --target benchmarktests_reduced_lung --parallel 4
```

Run only the distributed reduced-lung benchmark entries with MPI:

```text
mpirun -np 2 ./build/debug/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/Distributed
```

Useful subsets:

```text
mpirun -np 2 ./build/debug/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/Distributed/FullSolve
mpirun -np 2 ./build/debug/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/Distributed/LinearSolve
```

## Current Limitations

- The new YAML file is registered as an NP2 CTest input.
- The distributed benchmark is added to the existing benchmark executable, but distributed benchmark execution still needs to be launched manually with `mpirun`.
- The benchmark physical model is limited to rigid linear airways.
- The benchmark does not yet include Kelvin-Voigt airways, nonlinear airway resistance, terminal units, or mixed airway/terminal-unit trees.
- The sparse distributed comparison still uses the configured solver parameters. In this benchmark fixture the selected solver is UMFPACK, matching the existing reduced-lung benchmark convention.
- Build, distributed benchmark, and manual NP2 runtime verification were run after adding the file.

## Intended Next Step

The next logical step is to add an NP4 runtime variant if needed, and then extend distributed benchmarks to mixed airway/terminal-unit cases.
