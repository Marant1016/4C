# Reduced Lung Residual Evaluation Optimization Guide

Date: 2026-08-06

## Scope

Optimize residual evaluation in the custom reduced-lung `NewtonTree` solver, focusing initially on large serial tree cases such as:

```text
reduced_lung_lung_tree_gen16_500steps_newton_tree.4C.yaml
```

## Baseline

Current gen16 profile:

```text
newton_total_s: 13.2191
residual_s: 5.94365
tree_assembly_s: 1.17429
tree_solve_s: 4.04021
```

Residual evaluation is currently the largest measured part of the Newton solve.

The active gen16 models are:

* Rigid airways
* Linear flow resistance
* No airway inertia
* Linear terminal-unit elasticity
* Kelvin-Voigt terminal-unit rheology
* Zero terminal-unit viscosity
* Pressure boundary conditions

## Goals

* Reduce `residual_s` without changing model equations or Newton convergence behavior.
* Remove unnecessary allocations, temporary vectors, repeated API calls, and copies.
* Optimize the common gen16 physics while preserving generic scalar fallback paths.
* Introduce SIMD only after simpler scalar optimizations have been measured.

## Step 1: Add Residual Phase Profiling

Split residual timing into:

```text
residual_clear_s
residual_airways_s
residual_terminal_units_s
residual_junctions_s
residual_boundary_conditions_s
residual_other_s
residual_norm_s
residual_evaluations
```

Add named residual assembler phases to the assembly pipeline and time:

* residual clearing;
* each model callback group;
* residual norm calculation.

This profiling must be completed before changing residual algorithms.

## Step 2: Use Cached Spans and Direct Local Writes

Replace per-entry calls such as:

```cpp
target.replace_local_value(row, value);
```

with direct writes to a cached local residual span where the row is local and uniquely owned.

Cache before each hot loop:

* DOF value span;
* residual value span;
* local row IDs;
* local DOF ID vectors;
* frequently used model arrays.

Apply this first to:

* airway wall residuals;
* terminal-unit rheology residuals;
* junction residuals;
* boundary-condition residuals.

Keep the existing API where row locality or ownership is uncertain.

## Step 3: Specialize the Common Airway Residual

Add a direct residual evaluator for:

```text
Rigid wall + linear resistance + no inertia
```

For this path:

* precompute constant Poiseuille resistance once;
* evaluate the complete residual in one loop;
* avoid resistance and inertia scratch vectors;
* write results directly into residual storage.

Keep the generic implementation for nonlinear resistance, inertia, Kelvin-Voigt walls, and unsupported model combinations.

## Step 4: Fuse the Common Terminal-Unit Residual

Add a fused evaluator for:

```text
Linear elasticity + Kelvin-Voigt rheology
```

Compute the elastic pressure and rheology residual in one loop instead of writing elastic pressure to a temporary vector and reading it again.

For the gen16 case, where `eta = 0`, use a dedicated path that avoids the viscosity calculation.

Preserve stored elastic-pressure values where they are required for output.

## Step 5: Remove Scratch Work From Hot Paths

Audit persistent and temporary residual vectors.

Rules:

* Precompute values that remain constant.
* Compute single-use dynamic values directly in the consumer loop.
* Retain persistent storage only when values are reused elsewhere.
* Avoid repeated `resize`, `clear`, fill, and copy operations.

For gen16:

* airway resistance is constant;
* inertia is always zero;
* terminal elastic pressure can be fused into residual evaluation.

## Step 6: Reprofile Before Adding SIMD

After the scalar optimizations, rerun profiling.

Add batch/SIMD kernels only for phases that remain significant. Initial candidates are:

* rigid linear airway residuals;
* linear Kelvin-Voigt terminal-unit residuals;
* bifurcation junction residuals;
* pressure boundary conditions.

Use SoA model arrays, padded or masked tails, and scalar fallback paths.

## Step 7: Revisit Clear and Norm Costs

Only optimize these operations if profiling shows that they have become significant:

```cpp
residual_.put_scalar(0.0);
residual_.norm_2(&residual_norm);
```

Remove residual clearing only after proving that every residual row is overwritten exactly once.

Use a direct local norm only for serial execution. Distributed execution must retain the correct global reduction.

Do not remove the final convergence residual evaluation.

## Validation

After every implementation step, run:

```bash
git diff --check
cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

Add comparison tests for every specialized path:

* optimized residual versus generic residual;
* scalar versus SIMD residual;
* odd-sized SIMD batches;
* unsupported-model fallback behavior.

Scalar refactors should produce identical or near-machine-precision results.

## Benchmark

Use the release gen16 profile:

```bash
FOUR_C_REDUCED_LUNG_TREE_PROFILE=1 ./build/release/4C \
  ../files/reduced_lung_lung_tree_gen16_500steps_newton_tree.4C.yaml \
  ../output/gen16_newton_tree_profile
```

Run at least three repetitions under similar machine conditions.

Compare:

```text
Calculation
newton_total_s
residual_s
residual phase timers
tree_assembly_s
tree_solve_s
```

Expected result:

* `residual_s` decreases;
* `newton_total_s` decreases;
* tree assembly and tree solve times remain approximately unchanged.

## Implementation Order

1. Add residual phase profiling.
2. Introduce cached spans and direct local writes.
3. Specialize the rigid linear no-inertia airway residual.
4. Fuse the linear Kelvin-Voigt terminal-unit residual.
5. Remove unnecessary scratch operations.
6. Measure again.
7. Add SIMD only to remaining dominant phases.
8. Optimize residual clear or norm only if profiling justifies it.

## Expected End State

Residual evaluation should have substantially less API overhead, temporary storage, and repeated computation.

It is acceptable for residual evaluation to remain one of the largest components if its cost is comparable to the tree solve and no obvious unnecessary work remains.
