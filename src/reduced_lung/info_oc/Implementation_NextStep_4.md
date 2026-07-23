# Implementation NextStep 4 - Benchmark And Profile

This step adds optional profiling hooks and a reduced-lung Google Benchmark harness for comparing `Nox`, `NewtonSparse`, and `NewtonTree` on serial reduced-lung inputs.

No nonlinear-solver policy was changed. `Nox` remains the default/reference path, `NewtonTree` remains serial-only, and all profiling pointers default to `nullptr` so normal runtime behavior remains unchanged.

## Profiling API

Added:

- `src/reduced_lung/src/4C_reduced_lung_solver_profile.hpp`

New profile structs:

- `NoxSolverProfile`
- `NewtonSolverProfile`
- `SparseNewtonLinearSolverProfile`
- `TreeNewtonLinearSolverProfile`

The profile structs collect phase timings and counters when explicitly attached by benchmark code.

## Instrumented Code

Updated:

- `src/reduced_lung/src/4C_reduced_lung_helpers.hpp`
- `src/reduced_lung/src/4C_reduced_lung_helpers.cpp`
- `src/reduced_lung/src/4C_reduced_lung_newton_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_newton_solver.cpp`
- `src/reduced_lung/src/4C_reduced_lung_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_linear_solver.cpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

Measured phases/counters include:

- full NOX solve time
- NOX residual and sparse-Jacobian callback counts
- NOX residual callback assembly time
- NOX sparse-Jacobian callback assembly time
- NOX sparse matrix completion time
- custom Newton full solve time
- custom Newton state sync time
- custom Newton residual assembly time
- custom Newton sparse Jacobian assembly time
- custom Newton sparse matrix completion time
- custom Newton structured tree-linearization assembly time
- custom Newton linear correction solve time
- sparse linear solver solve time
- tree solver total solve time
- tree bottom-up condensation time
- tree top-down recovery time
- tree dense local solve time and count
- tree coefficient lookup time and count
- tree symbolic workspace counters: element count, total local block dofs, max local block size
- last custom Newton residual norm history and increment norm history

## Benchmark Harness

Added:

- `src/reduced_lung/benchmark_tests/CMakeLists.txt`
- `src/reduced_lung/benchmark_tests/4C_reduced_lung_solver_benchmark.cpp`

The benchmark target is registered through the existing project helper:

```text
four_c_auto_define_benchmark_tests()
```

Benchmark groups:

- `ReducedLung/FullSolve/...`: full nonlinear solve comparison for `Nox`, `NewtonSparse`, and `NewtonTree`.
- `ReducedLung/AssemblyPhases/BalancedAirways`: residual assembly, sparse Jacobian assembly, sparse matrix completion, and structured tree-linearization assembly measured separately.
- `ReducedLung/LinearSolve/BalancedAirways/Sparse`: sparse linear-solve cost for growing balanced airway trees.
- `ReducedLung/LinearSolve/BalancedAirways/StructuredTree`: structured tree-solve cost, local dense solve cost, coefficient lookup cost, and workspace counters for growing balanced airway trees.

Benchmark cases:

- single terminal unit
- 9-element serial airway chain
- 4-level balanced airway tree
- balanced airway trees with 2, 3, 4, and 5 levels for linear-solve scaling

## Benchmark Run

The benchmark target was built with the CLion-bundled CMake 4.2.2 because the shell-default `/usr/bin/cmake` is CMake 3.28.3 and the project requires CMake >= 3.30.

Build command used:

```text
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target benchmarktests_reduced_lung --parallel 4
```

Direct benchmark validation command used:

```text
./build/debug/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung --benchmark_min_time=0.01s
```

The direct run completed all 18 reduced-lung benchmark entries. Selected reported entries from the debug build:

```text
ReducedLung/FullSolve/SingleTerminalUnit/Nox                  15828 us
ReducedLung/FullSolve/SingleTerminalUnit/NewtonSparse          6803 us
ReducedLung/FullSolve/SingleTerminalUnit/NewtonTree             147 us
ReducedLung/FullSolve/SerialAirways/Nox                       11572 us
ReducedLung/FullSolve/SerialAirways/NewtonSparse              10443 us
ReducedLung/FullSolve/SerialAirways/NewtonTree                  410 us
ReducedLung/FullSolve/BalancedAirways/Nox                     11395 us
ReducedLung/FullSolve/BalancedAirways/NewtonSparse            12134 us
ReducedLung/FullSolve/BalancedAirways/NewtonTree                597 us
ReducedLung/AssemblyPhases/BalancedAirways                      141 us
ReducedLung/LinearSolve/BalancedAirways/Sparse/2               3623 us
ReducedLung/LinearSolve/BalancedAirways/Sparse/3               4043 us
ReducedLung/LinearSolve/BalancedAirways/Sparse/4               5927 us
ReducedLung/LinearSolve/BalancedAirways/Sparse/5               9000 us
ReducedLung/LinearSolve/BalancedAirways/StructuredTree/2       11.5 us
ReducedLung/LinearSolve/BalancedAirways/StructuredTree/3       21.7 us
ReducedLung/LinearSolve/BalancedAirways/StructuredTree/4        109 us
ReducedLung/LinearSolve/BalancedAirways/StructuredTree/5        118 us
```

The run emitted expected warnings about CPU scaling and debug-build timing noise. NOX also prints its existing solver status output during the full-solve benchmark cases.

## Allocation Profiling

The benchmark records tree workspace and operation counters directly:

- `tree_workspace_dofs`
- `tree_max_block`
- `tree_dense_solves`
- `tree_lookups`

External heap allocation profiling was not run in this step. Recommended commands after building the benchmark executable:

```text
heaptrack ./build/debug/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/LinearSolve/BalancedAirways/StructuredTree
valgrind --tool=massif ./build/debug/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/LinearSolve/BalancedAirways/StructuredTree --benchmark_min_time=0.01s
```

## Verification Run

Commands run successfully:

```text
cmake --build build/debug --target unittests_reduced_lung --parallel 4
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target benchmarktests_reduced_lung --parallel 4
./build/debug/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung --benchmark_min_time=0.01s
ctest -R "^unittests_reduced_lung$" --output-on-failure
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target unittests_reduced_lung.np2 --parallel 4
ctest -R "^unittests_reduced_lung(\.np2)?$" --output-on-failure
ctest -R "reduced_lung_.*newton_.*\.4C\.yaml-p1$" --output-on-failure
git diff --check
```

The runtime YAML regex also selected `test_cleanup`; all selected tests passed.

## How To Run The Benchmark

From the repo root:

```text
cd /scratch/Rodriguez/workspace/4C/4C
```

Use the CLion CMake, not `/usr/bin/cmake`, because `/usr/bin/cmake` is version 3.28.3 on this machine:

```text
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake \
  -S /scratch/Rodriguez/workspace/4C/4C \
  -B /scratch/Rodriguez/workspace/4C/4C/build/debug \
  -DFOUR_C_WITH_GOOGLE_BENCHMARK=ON \
  -DFOUR_C_ENABLE_FULL_BENCHMARK_TESTS=ON
```

Build only the reduced-lung benchmark executable:

```text
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake \
  --build /scratch/Rodriguez/workspace/4C/4C/build/debug \
  --target benchmarktests_reduced_lung \
  --parallel 4
```

Run all reduced-lung benchmarks directly:

```text
./build/debug/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung
```

For a quick validation run in debug mode:

```text
./build/debug/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung --benchmark_min_time=0.01s
```

Useful subsets:

```text
./build/debug/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/FullSolve --benchmark_min_time=0.01s
./build/debug/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/AssemblyPhases --benchmark_min_time=0.01s
./build/debug/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/LinearSolve --benchmark_min_time=0.01s
```

CTest can also run the benchmark target:

```text
ctest --test-dir build/debug -R "^benchmarktests_reduced_lung$" --output-on-failure
```

## Remaining After This Step

- Benchmark timings shown above are from a debug build and should not be used as final performance numbers.
- A release build should be used for thesis-quality measurements.
- External heap profiling still needs to be run separately.
- Parallel tree solving remains future work.
