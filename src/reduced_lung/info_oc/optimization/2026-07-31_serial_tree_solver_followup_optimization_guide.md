# Serial Tree Solver Follow-Up Optimization Guide

Date: 2026-07-31

## Scope

This guide defines the next reduced-lung serial `TreeNewtonLinearSolver` optimization steps after the Step 2-8 grouped/vectorization-friendly refactor.

Only the serial tree solver path is in scope unless a later step explicitly says otherwise. The distributed tree solver should remain unchanged while these serial-only follow-up optimizations are evaluated.

## Current Benchmark Signal

The Step 8 rerun showed mixed results:

- Standalone `StructuredTree/5_mean` linear solve was essentially flat versus the Step 1 baseline.
- Larger structured linear cases showed much lower `tree_dense_s`, which confirms that the batch dense-solver work is visible.
- `tree_lookup_s` remains a large cost and increased in the standalone structured linear-solve rerun.
- `tree_top_down_s` increased in the standalone structured linear-solve rerun, especially for small trees where grouping/staging overhead dominates.
- Full `BalancedAirways/NewtonTree_mean` improved substantially, but the CPU-scaling warning and benchmark CV mean it should be treated as directional, not final proof.

The follow-up work should therefore focus on lookup overhead and top-down staging overhead before adding more broad restructuring.

## Measurement Rule For All Follow-Up Steps

Before and after each optimization step, run longer release benchmarks than the Step 1/Step 8 quick runs when feasible:

```text
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/release --target benchmarktests_reduced_lung --parallel 4
./build/release/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/LinearSolve/BalancedAirways/StructuredTree --benchmark_repetitions=5 --benchmark_min_time=0.1s
./build/release/tests/benchmarktests_reduced_lung --benchmark_filter="ReducedLung/FullSolve/(SingleTerminalUnit|SerialAirways|BalancedAirways)/NewtonTree" --benchmark_repetitions=5 --benchmark_min_time=0.1s
```

If possible, also pin the benchmark process to a stable CPU and disable CPU scaling. If that is not possible, document the warning and treat small changes as inconclusive.

## Step 9 - Precompute Coefficient Locations

### Objective

Remove repeated row-entry searches from the serial tree solver hot path.

### Reason

`tree_lookup_s` remains one of the largest measured costs. The current structured coefficient access still spends time finding values that are requested repeatedly with the same symbolic pattern.

### Implementation Direction

- During symbolic setup, identify every matrix coefficient needed by each element's bottom-up and top-down work.
- Store direct coefficient locations or compact access descriptors in the serial plan.
- During numeric solves, load coefficients through those descriptors instead of searching rows by global dof each time.
- Keep sparse and structured coefficient sources behaviorally equivalent.

### Validation

- Run the reduced-lung unit tests.
- Run the `newton_tree` reduced-lung input tests.
- Compare Step 9 benchmark counters against Step 8, especially `tree_lookup_s`, `tree_solve_s`, and full-solve `newton_total_s`.

### Success Criteria

- `tree_lookup_s` decreases for `StructuredTree/4_mean`, `StructuredTree/5_mean`, and `BalancedAirways/NewtonTree_mean`.
- No increase in nonlinear iteration counts.
- No change in final residual behavior beyond normal floating-point noise.

## Step 10 - Direct Tree-Block Coefficient Storage

### Objective

Store the coefficients needed by the serial tree solve in compact per-element tree blocks before solving.

### Reason

Precomputed locations still read through the original matrix/coefficient representation. Direct tree-block storage should improve locality and make bottom-up group assembly less dependent on global matrix layout.

### Implementation Direction

- During tree assembly or before the numeric solve, populate compact coefficient arrays per element.
- Store values in the order consumed by grouped bottom-up assembly.
- Keep the existing global matrix behavior intact for code outside the serial tree solver.
- Prefer a minimal representation first; do not introduce a new general matrix abstraction unless profiling proves it is needed.

### Validation

- Same correctness tests as Step 9.
- Compare `tree_lookup_s`, `tree_bottom_up_s`, and `tree_assembly_s` separately, because this step may shift cost from solve to assembly.

### Success Criteria

- Net `linear_solve_s` or full `newton_total_s` improves, not just `tree_lookup_s`.
- Any extra `tree_assembly_s` is smaller than the solve-time savings.

## Step 11 - Reduce Top-Down Staging Overhead

### Objective

Recover the top-down regressions introduced by grouped lane staging in Step 8.

### Reason

`tree_top_down_s` increased in the standalone structured linear-solve rerun. The grouped top-down path is more regular, but for tiny groups the extra gathers, buffers, and scatters can cost more than the scalar work it replaced.

### Implementation Direction

- Add a scalar or low-overhead fast path for single-lane and very small top-down groups.
- Select scalar versus grouped recovery during symbolic setup or at group execution time using a simple threshold.
- Keep the current grouped path for larger same-shape groups.
- Avoid adding backward-compatibility layers; this is an internal solver policy decision.

### Validation

- Compare `tree_top_down_s` for all structured tree sizes.
- Confirm `StructuredTree/2_mean` and `StructuredTree/3_mean` do not regress further.
- Confirm larger cases keep any benefit from grouped traversal.

### Success Criteria

- `tree_top_down_s` drops versus Step 8 for small and medium structured-tree cases.
- `tree_solve_s` does not regress for `StructuredTree/5_mean` or `BalancedAirways/NewtonTree_mean`.

## Step 12 - Specialize Common Element Shapes

### Objective

Reduce branch and loop overhead in the hottest repeated element shapes.

### Reason

The serial tree solver mostly sees a small set of local shapes: leaf elements, internal bifurcation elements, and root-like cases. Dedicated kernels for these shapes can remove repeated conditionals and make compiler optimization easier.

### Implementation Direction

- Profile or inspect group metadata to identify the most frequent shapes.
- Add dedicated kernels only for shapes that dominate the benchmark cases.
- Keep generic fallback logic for uncommon shapes.
- Start with the lowest-risk cases, such as leaf `2x2` and common internal two-child groups.

### Validation

- Same correctness tests as earlier steps.
- Compare `tree_bottom_up_s`, `tree_dense_s`, and `tree_top_down_s` by benchmark size.

### Success Criteria

- Larger grouped cases improve without making small cases worse.
- Code duplication remains contained and easy to review.

## Step 13 - Select Scalar Versus Grouped Solver By Tree Size

### Objective

Use the grouped/vectorization-friendly solver only when the tree is large enough to amortize its overhead.

### Reason

The Step 8 data suggests that small trees can regress because the grouped path performs more staging work than the old scalar path. A simple internal policy can preserve small-case performance while keeping the grouped path for larger trees.

### Implementation Direction

- During symbolic setup, compute cheap metrics such as element count, maximum group size, and number of non-leaf groups.
- Choose a scalar fast path for very small trees and the grouped path for larger trees.
- Keep the choice internal to `TreeNewtonLinearSolver` and record enough counters to understand which path was used.

### Validation

- Confirm `SingleTerminalUnit/NewtonTree_mean` and `StructuredTree/2_mean` improve versus Step 8.
- Confirm larger cases do not lose the grouped path's dense-solver and traversal benefits.

### Success Criteria

- Small-tree benchmarks recover toward Step 1 baseline behavior.
- Large-tree benchmarks remain flat or improve versus Step 8.

## Recommended Order

1. Step 9: Precompute coefficient locations.
2. Step 11: Reduce top-down staging overhead.
3. Step 10: Direct tree-block coefficient storage, if Step 9 shows lookup remains significant.
4. Step 13: Scalar versus grouped solver selection, if small-case regressions remain after Step 11.
5. Step 12: Shape specialization, only after benchmark data identifies specific shapes worth specializing.

## Documentation Expectations

For each implemented follow-up step, create a matching `Optimization_Step_N.md` note with:

- Scope and files changed.
- Exact behavior-preservation statement.
- Verification commands and results.
- Release benchmark command and mean rows.
- Delta versus the previous step and versus Step 1 when useful.
- A short conclusion explaining whether the change should be kept, revised, or followed by another optimization.
