# Serial Tree Newton Full-Batch SIMD Guide

Date: 2026-08-04

## Goal

Finish the serial structured `NewtonTree` workflow so large trees are processed through a full batch/SIMD path for all production element shapes, instead of falling back to per-element scalar work in the hot bottom-up and top-down solver phases.

The target solver is the serial `TreeNewtonLinearSolver` selected when `comm_size == 1` and `nonlinear_solver = NewtonTree`.

Primary files:

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp
```

Supporting files:

```text
src/reduced_lung/src/4C_reduced_lung_solver_profile.hpp
src/reduced_lung/src/4C_reduced_lung_tree_linearization.cpp
src/reduced_lung/src/4C_reduced_lung_tree_linearization.hpp
src/reduced_lung/src/CMakeLists.txt
src/reduced_lung/tests/4C_reduced_lung_tree_linear_solver_test.cpp
```

The distributed tree solver is out of scope for this guide.

## Current State

The serial tree solver already has substantial SIMD/batch infrastructure:

- `tree_solver_simd` wraps `std::experimental::native_simd<double>` when `<experimental/simd>` is available.
- `build_symbolic_plan()` groups elements by tree layer, local block size, and child count.
- The large-tree bottom-up path dispatches grouped work through `bottom_up_layer_groups_`.
- The large-tree top-down path dispatches grouped work through `top_down_layer_groups_`.
- `solve_2x2_batch(...)` and `solve_3x3_batch(...)` already perform SIMD arithmetic for full chunks.
- Structured coefficient values are cached in direct arrays after `set_tree_linearization(...)`.

The current solver is not yet a full batch workflow because several hot paths still loop element-by-element or use scalar tails/fallbacks.

## Production Cases To Cover

Current reduced-lung element layouts produce these block sizes:

- `block_size == 2`: rigid airways and terminal units with `num_dofs == 3`.
- `block_size == 3`: compliant/Kelvin-Voigt airways with `num_dofs == 4`.

Current tree topology supports these child counts:

- `child_count == 0`: leaf element with outlet boundary.
- `child_count == 1`: serial connection.
- `child_count == 2`: bifurcation.

Full serial batch coverage means explicit batch/SIMD paths for all six cases:

```text
2x2 + 0 child
2x2 + 1 child
2x2 + 2 children
3x3 + 0 child
3x3 + 1 child
3x3 + 2 children
```

Generic scalar support for future `block_size > 3` may remain, but it should not be used by current production cases.

## Implementation Principles

1. Keep the dependency order by tree layer.
2. Vectorize across independent elements inside one `ElementGroup`.
3. Restrict the optimized path to `TreeNewtonLinearSolverCoefficientSource::StructuredTreeBlocks`.
4. Keep scalar code for validation, small trees, unsupported future block sizes, and hard numerical fallbacks.
5. Add profile counters so performance runs can prove which path was used.

Do not try to vectorize the sparse-Jacobian coefficient source first. The sparse provider uses indirect row lookups and should stay a validation path unless there is a specific need to optimize it.

## Step 1: Add Explicit Batch Path Accounting

Add counters to `TreeNewtonLinearSolverProfile` in `4C_reduced_lung_solver_profile.hpp`.

Suggested counters:

```cpp
std::uint64_t simd_group_count = 0;
std::uint64_t simd_lane_count = 0;
std::uint64_t scalar_group_count = 0;
std::uint64_t scalar_tail_lane_count = 0;
std::uint64_t dense_fallback_count = 0;
std::uint64_t unsupported_block_fallback_count = 0;
```

Use these in the serial solver to distinguish:

- Full SIMD chunks.
- Scalar tails.
- Dense fallback solves.
- Unsupported shape fallbacks.

Acceptance criteria:

- Existing tests still pass.
- Profiling can report whether large-tree solves are actually using the batch path.

## Step 2: Replace Scalar Tails With Masked Or Padded SIMD

Current SIMD kernels process full chunks and then scalar tails. Replace tail handling in the hot paths with one of these approaches:

- Masked SIMD chunks for the final partial group.
- Padded group storage so each group size rounds up to `tree_solver_simd::width()`.

The padded approach is simpler because current `std::experimental::simd` support is implementation-dependent for masked stores.

Required changes:

- Add padded group size helpers in `4C_reduced_lung_tree_linear_solver.cpp`.
- Allocate batch workspaces using padded group sizes in `build_symbolic_plan()`.
- Store a lane-valid mask or valid count for the final chunk.
- Ensure padded lanes never write correction values, inlet pressure stamps, child pressures, or subtree relations.

Current scalar tails to remove or isolate:

```text
solve_2x2_batch(...)
solve_3x3_batch(...)
assemble_2x2_leaf_group(...)
assemble_2x2_one_child_group(...)
assemble_2x2_two_child_group(...)
recover_2x2_top_down_group(...)
```

Acceptance criteria:

- Large groups no longer increment `scalar_tail_lane_count` in normal runs.
- Numerical results match the sparse Newton solver and existing tree tests.

## Step 3: Convert 2x2 Batch Workspaces To Direct SoA

The existing `2x2` path uses SIMD arithmetic but still gathers from interleaved per-element workspace arrays. Add explicit SoA arrays for `2x2`, similar to existing `3x3` arrays.

Add to `TreeNewtonLinearSolver` in `4C_reduced_lung_tree_linear_solver.hpp`:

```cpp
std::vector<double> batch_2x2_a00_;
std::vector<double> batch_2x2_a01_;
std::vector<double> batch_2x2_a10_;
std::vector<double> batch_2x2_a11_;
std::vector<double> batch_2x2_rhs_constant0_;
std::vector<double> batch_2x2_rhs_constant1_;
std::vector<double> batch_2x2_rhs_inlet_pressure0_;
std::vector<double> batch_2x2_rhs_inlet_pressure1_;
std::vector<double> batch_2x2_intercept0_;
std::vector<double> batch_2x2_intercept1_;
std::vector<double> batch_2x2_slope0_;
std::vector<double> batch_2x2_slope1_;
```

Then update `solve_2x2_batch(...)` to operate on SoA arrays directly.

Acceptance criteria:

- `2x2` assembly can write directly into SoA arrays.
- `write_subtree_relation_group(...)` can read `batch_2x2_slope*` and `batch_2x2_intercept*` without reading interleaved `workspace_*` arrays.

## Step 4: Add 3x3 Structured Batch Assembly

The biggest current gap is that `3x3` bottom-up groups still use generic scalar `assemble_group(...)`.

Add these functions or lambdas near the existing `2x2` assembly group code:

```text
assemble_3x3_leaf_group(...)
assemble_3x3_one_child_group(...)
assemble_3x3_two_child_group(...)
```

Each function should mirror the existing `2x2` structured chunk pattern:

- Load three equation RHS values per lane.
- Load three inlet-pressure coefficients per lane.
- Load nine matrix coefficients per lane.
- For child slots, load pressure parent/child coefficients, pressure RHS, child subtree relation, and child-flow coefficient.
- Compute child pressure relation and downstream flow contribution in SIMD.
- Store child pressure relation and update the downstream flow equation row.

Only use the SIMD structured path when:

```cpp
coefficient_source_ == TreeNewtonLinearSolverCoefficientSource::StructuredTreeBlocks
```

The sparse coefficient source may keep the scalar generic `assemble_group(...)` path.

Acceptance criteria:

- `3x3 + 0 child`, `3x3 + 1 child`, and `3x3 + 2 children` no longer call generic scalar `assemble_group(...)` in structured large-tree runs.
- Kelvin-Voigt airway tree tests still match sparse Newton results.

## Step 5: Make 3x3 Solve Fully Batch-Oriented

`solve_3x3_batch(...)` already has SIMD arithmetic, but it first copies interleaved workspace data into SoA arrays with a scalar loop.

After Step 4, route the `3x3` structured assembly directly into the existing SoA arrays:

```text
batch_3x3_a00_ ... batch_3x3_a22_
batch_3x3_rhs_constant0_ ... batch_3x3_rhs_constant2_
batch_3x3_rhs_inlet_pressure0_ ... batch_3x3_rhs_inlet_pressure2_
```

Then update `solve_3x3_batch(...)` so the structured path skips the scalar copy loop.

Keep a scalar copy path only for sparse-Jacobian validation or unsupported block layouts.

Acceptance criteria:

- Structured `3x3` groups do not copy from `workspace_matrix_` into SoA arrays lane-by-lane before solving.
- The existing fallback path remains available for numerical failures.

## Step 6: Add 3x3 Top-Down SIMD Recovery

Current `recover_2x2_top_down_group(...)` has a SIMD chunk path. Current generic `recover_top_down_group(...)` is scalar lane-loop based and handles `3x3` groups.

Add:

```text
recover_3x3_top_down_group(...)
```

It should process all lanes in SIMD chunks:

- Load inlet pressure for each element.
- Compute three unknown correction values.
- Write inlet pressure correction and unknown corrections for valid lanes.
- Compute outlet pressure using `outlet_pressure_unknown_index_`.
- For each child slot, compute child inlet pressure and call `set_inlet_pressure(...)` for valid lanes.

Dispatch in top-down traversal:

```cpp
if (group.block_size == 2 && group.child_count <= 2)
{
  recover_2x2_top_down_group(group);
}
else if (group.block_size == 3 && group.child_count <= 2)
{
  recover_3x3_top_down_group(group);
}
else
{
  recover_top_down_group(group);
}
```

Acceptance criteria:

- Compliant airway groups no longer use scalar generic top-down recovery in large-tree structured runs.
- Child inlet pressures are stamped exactly once, as before.

## Step 7: Batch Subtree Relation Writes

Current `write_subtree_relation_group(...)` loops over grouped elements and reads interleaved workspace arrays.

Add specialized writers:

```text
write_2x2_subtree_relation_group(...)
write_3x3_subtree_relation_group(...)
```

These should use the SoA batch solve outputs and `inlet_flow_unknown_index_` to write:

```text
subtree_relation_G_[element] = slope[inlet_flow_unknown_index]
subtree_relation_h_[element] = intercept[inlet_flow_unknown_index]
```

This step is less arithmetic-heavy than solves, but it removes another per-element loop from the normal structured batch path.

Acceptance criteria:

- Normal `2x2` and `3x3` structured large-tree groups use specialized subtree relation writers.

## Step 8: Gate Or Remove Small-Group Scalar Bypasses

The solver intentionally uses scalar traversal for small trees and tiny top-down groups:

```text
use_scalar_tree_solve_ = element_count <= 7
```

For a big-tree-first workflow, these can remain. For strict batch behavior, make them configurable or disable them when a new profiling/debug option requests forced batch execution.

Recommended approach:

- Keep scalar small-tree behavior for production default.
- Add a developer-only forced-batch flag if needed for testing.
- Ensure big-tree benchmarks do not hit scalar thresholds.

Acceptance criteria:

- Existing small-tree tests remain stable.
- Big-tree tests and benchmarks show zero small-tree scalar traversal.

## Step 9: Define Numerical Fallback Policy

The existing batched dense solves fall back to scalar `solve_dense_system(...)` for unsafe pivots or non-finite values.

Policy options:

- Keep fallback as exceptional correctness behavior and count it.
- Implement batched pivoting for `3x3` so fallback is only for truly singular systems.

Recommended first implementation:

- Keep scalar fallback.
- Add `dense_fallback_count`.
- Treat nonzero fallback count in standard big-tree benchmark cases as a performance bug to investigate.

Acceptance criteria:

- Fallback count is visible in profiles.
- Normal validation inputs do not rely on scalar dense fallback.

## Step 10: Update Tests

Add or extend tests in:

```text
src/reduced_lung/tests/4C_reduced_lung_tree_linear_solver_test.cpp
```

Required test coverage:

- Big rigid airway tree covering `2x2 + 0/1/2 child`.
- Big Kelvin-Voigt airway tree covering `3x3 + 0/1/2 child`.
- Mixed airway/terminal-unit tree covering mixed `2x2` and `3x3` groups.
- Group sizes larger than SIMD width.
- Group sizes not divisible by SIMD width.
- Reused structured tree solver across at least two solves.

Each test should compare against sparse Newton or the existing sparse linear solver result within the existing tolerances.

Optional profile assertions:

- `simd_group_count > 0`.
- `simd_lane_count > 0`.
- `unsupported_block_fallback_count == 0` for current production block sizes.
- `dense_fallback_count == 0` for well-conditioned validation cases.

## Verification Commands

Use the project’s normal build/test commands for reduced lung. If a local build directory exists, typical verification targets are:

```text
git diff --check
cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

Also run at least one larger benchmark-style input, such as the generation-10 Newton tree input, and inspect profile counters to confirm the batch path is used.

## Expected End State

For large serial structured `NewtonTree` runs with current production element models:

- Bottom-up assembly is batched for `2x2` and `3x3` groups.
- Dense local solves are batched for `2x2` and `3x3` groups.
- Subtree relation writes are batched for `2x2` and `3x3` groups.
- Top-down recovery is batched for `2x2` and `3x3` groups.
- Scalar per-element loops remain only for small-tree mode, unsupported future block sizes, sparse-Jacobian validation, or exceptional numerical fallback.
