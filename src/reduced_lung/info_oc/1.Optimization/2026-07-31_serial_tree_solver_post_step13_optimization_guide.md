# Serial Tree Solver Post-Step-13 Optimization Guide

Date: 2026-07-31

## Scope

This guide defines the next serial `TreeNewtonLinearSolver` optimization sequence after `Optimization_Step_13.md`.

The focus is the new tree-based serial solver and its structured-tree workflow. The main goal is to make the batch-oriented path as fast as possible while preserving the current behavior and keeping `NewtonTree` an explicit opt-in solver.

Unless a future step explicitly says otherwise, leave the distributed tree solver unchanged.

## Current Benchmark Signal

Step 13 removed structured coefficient lookup cost from the serial solve loop and improved the standalone tree solve substantially. The remaining measured costs point to three major areas:

- Structured tree linearization assembly and coefficient refresh now matter at full-solve level.
- Top-down recovery remains a large part of the serial tree solve.
- The grouped/batch path still performs staging and copy work before doing small dense solves.

Important Step 13 reference values:

| Benchmark | tree_solve_s | tree_bottom_up_s | tree_top_down_s | tree_dense_s | tree_assembly_s |
| --- | ---: | ---: | ---: | ---: | ---: |
| `StructuredTree/5_mean` | `5.93419 us` | `2.92978 us` | `2.86843 us` | `0.524934 us` | n/a |
| `BalancedAirways/NewtonTree_mean` | `8.89066 us` | `4.31666 us` | `4.00306 us` | `0.976806 us` | `9.62805 us` |

The `BalancedAirways/NewtonTree_mean` full solve shows `tree_assembly_s` slightly above `tree_solve_s`, so future work should optimize both structured assembly/refresh and the linear solve itself.

## Measurement Rule

Before and after each optimization step, run the same correctness checks used in Step 13:

```text
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
git diff --check
```

For performance comparisons, use release benchmarks with a longer minimum time when feasible:

```text
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/release --target benchmarktests_reduced_lung --parallel 4
./build/release/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/LinearSolve/BalancedAirways/StructuredTree --benchmark_repetitions=5 --benchmark_min_time=0.1s
./build/release/tests/benchmarktests_reduced_lung --benchmark_filter="ReducedLung/FullSolve/(SingleTerminalUnit|SerialAirways|BalancedAirways)/NewtonTree" --benchmark_repetitions=5 --benchmark_min_time=0.1s
```

If CPU scaling or display authorization warnings appear, document them and treat small benchmark deltas as inconclusive.

## Step 14 - Reuse Structured Linearization And Coefficient Locations

### Objective

Reduce full-solve `tree_assembly_s` by avoiding repeated allocation and repeated row-entry scans when the structured linearization dimensions and sparsity pattern are unchanged.

### Reason

The full Newton path currently rebuilds `TreeLinearization` every correction and then calls `set_tree_linearization()`, which refreshes direct coefficient buffers by scanning rows again.

Relevant code:

- `src/reduced_lung/src/4C_reduced_lung_newton_solver.cpp`: `assemble_tree_linearization_for_current_state()`
- `src/reduced_lung/src/4C_reduced_lung_tree_linearization.cpp`: `TreeLinearization::reset()`, `clear_values()`, `set_value()`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`: `resolve_structured_coefficient_locations()`

### Implementation Direction

- Use `TreeLinearization::clear_values()` instead of `reset(...)` when row and dof dimensions are unchanged.
- Preserve row vector capacity across Newton corrections.
- In `resolve_structured_coefficient_locations()`, first try the existing `structured_entry_index` and only scan the row if the index is invalid or no longer matches the expected dof.
- Keep the existing assertions that validate row/dof bounds and descriptor consistency.
- Keep sparse-Jacobian behavior unchanged.

### Validation Focus

- `tree_assembly_s` in full `NewtonTree` benchmarks.
- `newton_total_s`, especially `BalancedAirways/NewtonTree_mean`.
- No change in `tree_lookups`, which should remain zero for structured serial solves.

### Success Criteria

- Full-solve `tree_assembly_s` decreases measurably.
- `tree_solve_s` does not regress.
- Nonlinear iteration counts remain unchanged.

## Step 15 - Write Delta Through Precomputed Local Indices

### Objective

Reduce top-down recovery cost by replacing per-dof global replacement calls with direct local vector writes in the serial solver.

### Reason

The serial solver already asserts that all correction dofs are locally available, but top-down recovery still writes each correction through `delta.replace_global_value(...)`.

Relevant code:

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`: `set_delta_value(...)`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`: `recover_top_down_element(...)` and `recover_top_down_group(...)`
- `src/core/linalg/src/sparse/4C_linalg_vector.hpp`: `get_values()` and `local_values_as_span()`

### Implementation Direction

- Add serial-only correction write indices to the symbolic plan.
- Precompute local correction ids for each element dof during `build_symbolic_plan()` or from metadata that is already local in serial.
- Validate once that each target correction dof is locally owned.
- Replace `delta.replace_global_value(global_dof_id, value)` with `delta.get_values()[local_index] = value` in the serial tree top-down path.
- Keep the distributed tree solver unchanged.

### Validation Focus

- `tree_top_down_s` in structured linear-solve benchmarks.
- `linear_solve_s` and `tree_solve_s` in full `NewtonTree` benchmarks.
- Correctness of every recovered correction dof.

### Success Criteria

- `tree_top_down_s` decreases for `StructuredTree/4_mean`, `StructuredTree/5_mean`, and `BalancedAirways/NewtonTree_mean`.
- No correctness regressions in unit or input tests.

## Step 16 - Specialize Common Top-Down Recovery Shapes

### Objective

Remove grouped top-down staging overhead for common element shapes.

### Reason

`recover_top_down_group()` stages inlet pressure, unknown values, outlet pressure, and child pressure through temporary lane buffers. That makes the loop regular, but it adds gather/compute/scatter phases that can dominate small per-element work.

Relevant code:

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`: `recover_top_down_group(...)`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`: `recover_top_down_element(...)`

### Implementation Direction

- Add direct top-down kernels for the common grouped shapes, starting with `block_size == 2` and `child_count == 0`, `1`, and `2`.
- Compute correction values directly from `workspace_slope_`, `workspace_intercept_`, and inlet pressure.
- Write directly to delta using the local-index path from Step 15.
- Keep the generic grouped path for uncommon shapes.
- Keep the scalar small-group fallback, but retune it only after this step is benchmarked.

### Validation Focus

- `tree_top_down_s` and total `tree_solve_s` for structured tree levels 4 and 5.
- Single-terminal and small serial-airway cases, to ensure extra dispatch does not hurt tiny trees.

### Success Criteria

- `tree_top_down_s` improves without increasing `tree_bottom_up_s` or `tree_dense_s`.
- Code duplication stays limited to the hot common shapes.

## Step 17 - Remove Redundant Batch Dense-Solver Staging

### Objective

Reduce bottom-up grouped solve overhead by avoiding workspace-to-batch copies before `2x2` and `3x3` solves.

### Reason

The current grouped bottom-up path writes matrix and RHS entries into generic workspace arrays, then the batch solvers copy those values into separate structure-of-arrays lane buffers before solving.

Relevant code:

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`: `solve_2x2_batch(...)`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`: `solve_3x3_batch(...)`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`: `assemble_2x2_leaf_group(...)`, `assemble_2x2_one_child_group(...)`, `assemble_2x2_two_child_group(...)`

### Implementation Direction

- For common `2x2` groups, compute intercept and slope directly from local scalar values or direct coefficient buffers.
- Avoid writing a full generic matrix and then copying it into `batch_2x2_*` arrays when the shape is known.
- Keep the generic `solve_dense_system(...)` fallback for uncommon block sizes and singular-pivot fallback cases.
- Evaluate `3x3` only after `2x2` specialization shows a clear win.

### Validation Focus

- `tree_bottom_up_s` and `tree_dense_s`.
- Larger structured linear-solve benchmarks where group sizes are large enough to amortize dispatch.

### Success Criteria

- `tree_bottom_up_s` decreases for `StructuredTree/5_mean` and `BalancedAirways/NewtonTree_mean`.
- Dense fallback behavior remains correct for near-singular local blocks.

## Step 18 - Gate Profiling Timers And Retune Thresholds

### Objective

Remove unnecessary timer overhead from production paths and retune scalar/grouped crossover thresholds after the structural optimizations above.

### Reason

The kernels now run in microseconds or sub-microseconds, so timer calls can distort detailed profile benchmarks when profiling is enabled. Also, local-write and staging changes will move the scalar/grouped crossover points.

Relevant code:

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`: solve, bottom-up, top-down, and dense-solve timing
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`: `use_scalar_tree_solve_`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`: `top_down_scalar_group_threshold`

### Implementation Direction

- Ensure expensive timing is only taken when `profile_ != nullptr`.
- Consider a coarse-profile mode for production-level timing and a detailed-profile mode for micro-kernel counters only if needed.
- Retune `scalar_tree_element_threshold` after Steps 14-17.
- Retune `top_down_scalar_group_threshold` after Step 16.
- Do not tune thresholds before reducing the top-down and assembly overheads, because the current crossover points will likely change.

### Validation Focus

- Structured linear-solve cases with and without profile counters if practical.
- Full `NewtonTree` benchmark total time, not just profile subcounters.

### Success Criteria

- Production solve path does not pay for detailed profiling when no profile is attached.
- Threshold changes improve or preserve all Step 13 benchmark cases.

## Recommended Order

1. Step 14: Reuse structured linearization and coefficient locations.
2. Step 15: Write delta through precomputed local indices.
3. Step 16: Specialize common top-down recovery shapes.
4. Step 17: Remove redundant batch dense-solver staging.
5. Step 18: Gate profiling timers and retune thresholds.

## Documentation Expectations

For each implemented step, create a matching `Optimization_Step_N.md` file with:

- Scope and files changed.
- Behavior-preservation statement.
- Verification commands and results.
- Release benchmark commands and mean rows.
- Delta versus Step 13 and the immediately previous step.
- A short conclusion explaining whether the change should be kept, revised, or followed by another optimization.
