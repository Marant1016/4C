# SIMD Step 6

Date: 2026-08-03

## Goal

Add an explicit SIMD fast path to the serial `TreeNewtonLinearSolver` `3x3` dense batch solve.

This step targets:

```text
solve_3x3_batch(...)
```

and extends the dense-solve SIMD work beyond the `2x2` path from Step 3.

## Files Changed

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
```

## Implementation

Updated `solve_3x3_batch(...)` to use grouped offset caches from Step 2:

```text
grouped_unknown_begin
grouped_matrix_begin
```

instead of recomputing offsets through:

```text
element_index -> unknown_offset[element_index]
element_index -> matrix_offset[element_index]
```

The `3x3` call site now passes:

```text
grouped_unknown_begin_
grouped_matrix_begin_
```

## SIMD Fast Path

The existing `3x3` batch solver already stages matrix and RHS values into SoA arrays:

```text
batch_3x3_a00_ ... batch_3x3_a22_
batch_3x3_rhs_constant0_ ... batch_3x3_rhs_constant2_
batch_3x3_rhs_inlet_pressure0_ ... batch_3x3_rhs_inlet_pressure2_
```

The SIMD path processes full chunks of:

```text
tree_solver_simd::width()
```

It uses the same unpivoted `3x3` elimination formulas as the previous scalar fast path:

```text
p0 = a00
l10 = a10 / p0
l20 = a20 / p0

p1 = a11 - l10 * a01
u12 = a12 - l10 * a02
u21 = a21 - l20 * a01
u22 = a22 - l20 * a02

l21 = u21 / p1
p2 = u22 - l21 * u12
```

Then it solves both RHS systems with SIMD arithmetic:

```text
constant RHS -> intercept0, intercept1, intercept2
inlet-pressure RHS -> slope0, slope1, slope2
```

Accepted SIMD results are scattered to both:

```text
batch_3x3_intercept*_ / batch_3x3_slope*_
workspace_intercept_ / workspace_slope_
```

## Fallback Behavior

The SIMD path uses a conservative full-chunk fallback rule.

If any lane in a SIMD chunk has an invalid pivot, nonfinite intermediate, or nonfinite result, the entire chunk is sent to the existing pivoted dense fallback:

```text
solve_dense_system(..., block_size = 3, ...)
```

Validated quantities include:

```text
p0, p1, p2
l10, l20, l21
u12, u21, u22
intercept0, intercept1, intercept2
slope0, slope1, slope2
```

The scalar direct path remains active for tails after full SIMD chunks.

Scalar fallback remains active for:

- small trees using `use_scalar_tree_solve_`,
- non-`3x3` groups,
- scalar tails after full SIMD chunks,
- invalid SIMD chunks,
- builds without `<experimental/simd>`,
- pivoted `solve_dense_system(...)` fallback.

## Behavior Preserved

- The mathematical formulas match the previous scalar `3x3` fast path.
- Unsafe `3x3` systems still use the existing pivoted dense fallback.
- Group traversal order is unchanged.
- `2x2` assembly, dense solve, and top-down recovery are unchanged.
- Top-down generic recovery is unchanged.
- The distributed tree solver is unchanged.

## Expected Benefit

This step targets `tree_dense_s` for cases with meaningful `3x3` local blocks.

The benchmark impact may be limited if the workload is mostly `2x2`. It should matter more for reduced-lung configurations that produce many `3x3` local systems.

## Verification

Verification target for this step:

```text
git diff --check
cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

Performance check target after a release build is available:

```text
./build/release/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/LinearSolve/BalancedAirways/StructuredTree
```

Record at least:

```text
tree_solve_s
tree_bottom_up_s
tree_dense_s
tree_top_down_s
```

## Follow-Up

After Step 6, explicit SIMD covers the main common serial grouped paths:

```text
2x2 bottom-up assembly
2x2 dense solve
2x2 top-down recovery
3x3 dense solve
```

The next decision should be benchmark-driven. If gather/scatter overhead dominates, consider packed per-group SoA buffers or finer fallback handling for invalid SIMD lanes.
