# Serial Tree Solver SIMD Implementation Guide

Date: 2026-08-03

## Goal

Implement an explicit SIMD workflow for the serial `TreeNewtonLinearSolver` so same-shape elements in the same tree layer are processed in vector batches instead of one scalar element at a time.

The first target is the serial structured `NewtonTree` path in:

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp
```

Do not change the distributed tree solver in the first SIMD milestone.

## Current State

The serial tree solver is already grouped, but not explicitly SIMD-vectorized.

Current behavior:

- Small trees use scalar full-tree traversal when `element_count <= scalar_tree_element_threshold`.
- Large trees are grouped by `block_size` and `child_count` inside each tree layer.
- Bottom-up traversal dispatches `bottom_up_layer_groups_`.
- Top-down traversal dispatches `top_down_layer_groups_`.
- Group functions still use scalar loops over `grouped_index` or `lane`.
- No explicit SIMD API, intrinsic, or SIMD pragma is currently used in the serial tree solver.

Important existing grouped paths:

```text
build_layer_groups(...)
assemble_2x2_leaf_group(...)
assemble_2x2_one_child_group(...)
assemble_2x2_two_child_group(...)
solve_2x2_batch(...)
solve_3x3_batch(...)
recover_2x2_top_down_group(...)
recover_top_down_group(...)
```

## Core SIMD Strategy

Tree dependencies prevent SIMD across parent-child levels. SIMD must happen across independent elements within the same tree layer.

Correct batching axis:

```text
same tree layer
+ same block_size
+ same child_count
= SIMD batch across independent elements
```

This matches the existing grouping model. The SIMD implementation should process each `ElementGroup` in chunks of the native vector width.

Example structure:

```cpp
using Vec = std::experimental::native_simd<double>;
constexpr int width = static_cast<int>(Vec::size());

int grouped_index = group.begin;
for (; grouped_index + width <= group.end; grouped_index += width)
{
  // Load width elements.
  // Compute width local systems.
  // Scatter width results.
}

for (; grouped_index < group.end; ++grouped_index)
{
  // Scalar tail.
}
```

The goal is not to remove all loops. The goal is that each full loop iteration handles `width` elements.

## SIMD Backend

Use C++20-compatible experimental SIMD:

```cpp
#include <experimental/simd>

namespace stdx = std::experimental;
using SimdDouble = stdx::native_simd<double>;
using SimdMask = typename SimdDouble::mask_type;
```

Rationale:

- The repository is C++20.
- Local GCC 13 and Clang 18 accept `<experimental/simd>` in syntax-only checks.
- This avoids hard-coding AVX/SSE intrinsics in the first implementation.
- The compiler can choose the native SIMD width for the active target.

Build note:

```text
FOUR_C_ENABLE_NATIVE_OPTIMIZATIONS=ON
```

should be enabled for performance benchmarks so the compiler can emit hardware-specific vector instructions such as AVX/AVX2/AVX-512 where available.

## Scope Rules

First SIMD milestone includes:

- Serial `TreeNewtonLinearSolver` only.
- Structured coefficient source first: `TreeNewtonLinearSolverCoefficientSource::StructuredTreeBlocks`.
- `2x2` grouped bottom-up solve first.
- Scalar fallback for tails, small groups, singular systems, and unsupported shapes.

First SIMD milestone excludes:

- Distributed tree solver changes.
- Sparse-Jacobian coefficient source SIMD.
- New mathematical formulation.
- Reordering equations or unknowns.
- Removing scalar fallback.

## Why Structured Coefficients First

The sparse-Jacobian provider performs row lookups and indirect searches. That makes clean SIMD loading difficult and can dominate the arithmetic.

The structured path already has direct coefficient buffers such as:

```text
equation_inlet_pressure_coefficient_values_
matrix_coefficient_values_
child_pressure_parent_coefficient_values_
child_pressure_child_coefficient_values_
child_flow_coefficient_values_
```

These direct arrays make SIMD gather or prepacked SoA loading practical.

## Data Layout Requirement

Current grouped traversal still often starts from `grouped_element_indices_` and then loads per-element offsets:

```cpp
const int element_index = grouped_element_indices_[grouped_index];
const int unknown_begin = unknown_offset_[element_index];
const int matrix_begin = matrix_offset_[element_index];
const int child_begin = child_interface_offset_[element_index];
```

This works, but it creates indirect scalar setup per lane. Before deeper SIMD work, add grouped metadata arrays parallel to `grouped_element_indices_`:

```cpp
std::vector<int> grouped_unknown_begin_;
std::vector<int> grouped_matrix_begin_;
std::vector<int> grouped_equation_begin_;
std::vector<int> grouped_child_begin_;
std::vector<int> grouped_inlet_pressure_local_dof_;
std::vector<int> grouped_outlet_pressure_unknown_index_;
std::vector<int> grouped_inlet_flow_unknown_index_;
```

Populate these in `build_symbolic_plan()` immediately after pushing an element into `grouped_element_indices_`.

This keeps the SIMD kernels simple and reduces repeated indirection in hot loops.

## Step 1: Add SIMD Helper Layer

Add a small helper section in the anonymous namespace of `4C_reduced_lung_tree_linear_solver.cpp`.

Suggested helpers:

```cpp
#include <experimental/simd>

namespace stdx = std::experimental;
using SimdDouble = stdx::native_simd<double>;
using SimdMask = typename SimdDouble::mask_type;

template <typename Load>
SimdDouble gather_simd(int grouped_begin, Load&& load)
{
  return SimdDouble([&](auto lane)
  {
    return load(grouped_begin + static_cast<int>(lane));
  });
}

template <typename Store>
void scatter_simd(const SimdDouble& values, int grouped_begin, Store&& store)
{
  for (int lane = 0; lane < static_cast<int>(SimdDouble::size()); ++lane)
  {
    store(grouped_begin + lane, values[static_cast<std::size_t>(lane)]);
  }
}
```

The first implementation can use scalar scatter. That is acceptable because most wins should come from vectorizing arithmetic-heavy solve and recovery expressions.

Acceptance criteria:

- Code compiles with GCC and Clang in C++20 mode.
- No solver behavior changes yet.
- No distributed solver changes.

## Step 2: Add Grouped Offset Caches

Extend `TreeNewtonLinearSolver` with grouped metadata arrays.

In the header:

```cpp
std::vector<int> grouped_unknown_begin_;
std::vector<int> grouped_matrix_begin_;
std::vector<int> grouped_equation_begin_;
std::vector<int> grouped_child_begin_;
```

In `build_layer_groups(...)`, when adding an element to `grouped_element_indices_`, also add its hot offsets:

```cpp
grouped_element_indices_.push_back(element_index);
grouped_unknown_begin_.push_back(unknown_offset_[element_index_size]);
grouped_matrix_begin_.push_back(matrix_offset_[element_index_size]);
grouped_equation_begin_.push_back(equation_offset_[element_index_size]);
grouped_child_begin_.push_back(child_interface_offset_[element_index_size]);
```

Validation should assert all grouped arrays have the same size as `grouped_element_indices_`.

Acceptance criteria:

- Existing tests pass.
- Existing scalar/grouped behavior is unchanged.
- No measurable regression in structured tree benchmark.

## Step 3: SIMD `2x2` Dense Batch Solve

Start with `solve_2x2_batch(...)` because it is compact and arithmetic-heavy.

Current scalar formula:

```text
det = a00 * a11 - a01 * a10

intercept0 = (rhs_constant0 * a11 - a01 * rhs_constant1) / det
intercept1 = (a00 * rhs_constant1 - rhs_constant0 * a10) / det

slope0 = (rhs_inlet_pressure0 * a11 - a01 * rhs_inlet_pressure1) / det
slope1 = (a00 * rhs_inlet_pressure1 - rhs_inlet_pressure0 * a10) / det
```

SIMD kernel outline:

```cpp
void solve_2x2_batch_simd(...)
{
  const int width = static_cast<int>(SimdDouble::size());
  int grouped_index = group_begin;

  for (; grouped_index + width <= group_end; grouped_index += width)
  {
    const SimdDouble a00 = gather_simd(grouped_index, load_a00);
    const SimdDouble a01 = gather_simd(grouped_index, load_a01);
    const SimdDouble a10 = gather_simd(grouped_index, load_a10);
    const SimdDouble a11 = gather_simd(grouped_index, load_a11);
    const SimdDouble rhs0 = gather_simd(grouped_index, load_rhs0);
    const SimdDouble rhs1 = gather_simd(grouped_index, load_rhs1);
    const SimdDouble pin0 = gather_simd(grouped_index, load_pin0);
    const SimdDouble pin1 = gather_simd(grouped_index, load_pin1);

    const SimdDouble det = a00 * a11 - a01 * a10;
    const auto valid = stdx::isfinite(det) && (stdx::abs(det) > SimdDouble(pivot_tolerance));

    if (stdx::all_of(valid))
    {
      const SimdDouble inv_det = SimdDouble(1.0) / det;
      const SimdDouble intercept0 = (rhs0 * a11 - a01 * rhs1) * inv_det;
      const SimdDouble intercept1 = (a00 * rhs1 - rhs0 * a10) * inv_det;
      const SimdDouble slope0 = (pin0 * a11 - a01 * pin1) * inv_det;
      const SimdDouble slope1 = (a00 * pin1 - pin0 * a10) * inv_det;

      scatter_simd(intercept0, grouped_index, store_intercept0);
      scatter_simd(intercept1, grouped_index, store_intercept1);
      scatter_simd(slope0, grouped_index, store_slope0);
      scatter_simd(slope1, grouped_index, store_slope1);
    }
    else
    {
      // Use existing scalar path per lane for this chunk.
    }
  }

  // Existing scalar tail.
}
```

Important fallback rule:

```text
If any lane in a SIMD chunk fails the determinant check, either fallback the full chunk or fallback only invalid lanes.
```

For the first implementation, fallback the full chunk. It is simpler and safer. Later, optimize to fallback only invalid lanes.

Acceptance criteria:

- Numerical results match existing tests.
- `tree_dense_s` improves or remains neutral for large structured `2x2` trees.
- `tree_solve_s` does not regress for small cases.

## Step 4: SIMD `2x2` Structured Assembly

After the `2x2` solve kernel is SIMD, move to assembly. Prioritize structured coefficient source only.

Targets:

```text
assemble_2x2_leaf_group(...)
assemble_2x2_one_child_group(...)
assemble_2x2_two_child_group(...)
```

For leaf groups, SIMD mainly loads residual and direct coefficient values into workspace arrays.

For one-child and two-child groups, vectorize the child-substitution arithmetic:

```text
child_pressure_slope = -pressure_parent_coeff / pressure_child_coeff
child_pressure_intercept = pressure_rhs / pressure_child_coeff

child_flow_slope = subtree_relation_G[child] * child_pressure_slope
child_flow_intercept = subtree_relation_G[child] * child_pressure_intercept + subtree_relation_h[child]

matrix_last_row_outlet += flow_child_coeff * child_flow_slope
rhs_shift += flow_child_coeff * child_flow_intercept
```

Keep scalar fallback for sparse coefficient source and unsupported group shapes.

Acceptance criteria:

- `tree_bottom_up_s` improves on large structured tree benchmark.
- Sparse-Jacobian tree solver path remains correct.
- Tests pass with assertions enabled.

## Step 5: SIMD `2x2` Top-Down Recovery

Target:

```text
recover_2x2_top_down_group(...)
```

Vectorize these expressions across the group:

```text
value0 = slope0 * inlet_pressure + intercept0
value1 = slope1 * inlet_pressure + intercept1
outlet_pressure = outlet_slope * inlet_pressure + outlet_intercept
child_pressure = child_pressure_slope * outlet_pressure + child_pressure_intercept
```

Writes to `delta` and child inlet-pressure slots may remain scalar scatter operations.

Acceptance criteria:

- `tree_top_down_s` improves on large structured tree benchmark.
- Correction vector exactly matches scalar reference within existing tolerances.
- No top-down double-write validation failures.

## Step 6: SIMD `3x3` Dense Batch Solve

Only do this after the `2x2` SIMD path is correct and benchmarked.

The existing `solve_3x3_batch(...)` already uses SoA temporary arrays:

```text
batch_3x3_a00_
batch_3x3_a01_
batch_3x3_a02_
...
batch_3x3_rhs_constant0_
batch_3x3_rhs_constant1_
batch_3x3_rhs_constant2_
...
```

This is a good SIMD target because loads are contiguous by lane after staging.

Use the same safety rule:

```text
fast SIMD path when all pivots are finite and above tolerance
scalar pivoted fallback otherwise
```

Acceptance criteria:

- `3x3` benchmark cases improve or remain neutral.
- Unsafe lanes still use the pivoted scalar fallback.
- Existing correctness tests pass.

## Step 7: Optional Group Packing (IGNORED)

If gather-heavy SIMD does not improve performance enough, add per-group packed SoA work buffers.

Example packed arrays for a `2x2` group:

```text
packed_a00[i]
packed_a01[i]
packed_a10[i]
packed_a11[i]
packed_rhs0[i]
packed_rhs1[i]
packed_pin0[i]
packed_pin1[i]
```

Then the SIMD solve can use contiguous `copy_from(...)` loads instead of gather constructors.

Tradeoff:

- Packing costs extra memory writes.
- Packed contiguous loads are often faster and easier for the compiler.

Only add this if benchmark data shows gather/indirection dominates.

## Step 8: Optional SIMD Counters (IGNORED)

Add profiling counters only after the first SIMD path works.

Potential fields in `TreeNewtonLinearSolverProfile`:

```cpp
std::uint64_t simd_2x2_chunks = 0;
std::uint64_t simd_2x2_lanes = 0;
std::uint64_t simd_2x2_fallback_chunks = 0;
std::uint64_t simd_2x2_fallback_lanes = 0;
std::uint64_t scalar_tail_lanes = 0;
```

These are useful for confirming that large cases actually enter the SIMD path.

Do not add counters before the implementation is stable; profiling code can add noise to microbenchmarks.

## Step 9: Sparse Path Decision (IGNORED)

The sparse-Jacobian path should stay scalar unless there is a concrete benchmark reason to optimize it.

Reason:

- It performs sparse row extraction and column search.
- Memory access is irregular.
- Structured `NewtonTree` is the production performance path.

If needed later, optimize sparse by pre-extracting the same direct coefficient buffers used by the structured path before entering grouped traversal.

## Step 10: Benchmark Protocol

Correctness checks after each step:

```text
git diff --check
cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

Performance checks after each SIMD milestone:

```text
./build/release/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/LinearSolve/BalancedAirways/StructuredTree
./build/release/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/FullSolve/.*/NewtonTree
```

Record at least:

```text
tree_solve_s
tree_bottom_up_s
tree_top_down_s
tree_dense_s
tree_lookup_s
tree_dense_solves
tree_workspace_dofs
tree_max_block
```

For thesis-quality timing:

- Use a release build.
- Enable native optimization if allowed: `FOUR_C_ENABLE_NATIVE_OPTIMIZATIONS=ON`.
- Run multiple repetitions.
- Note CPU scaling warnings.
- Compare against the previous scalar/grouped implementation.

## Expected Risks

Risk: SIMD gather overhead cancels arithmetic gains.

Mitigation: add packed SoA buffers only after measuring.

Risk: invalid pivot handling diverges from scalar solver.

Mitigation: keep the existing pivoted scalar fallback and use it for unsafe chunks.

Risk: small trees regress due to SIMD setup overhead.

Mitigation: keep `use_scalar_tree_solve_` and scalar tails.

Risk: compiler support for `<experimental/simd>` differs across systems.

Mitigation: keep SIMD implementation isolated and add a compile-time fallback to scalar kernels if needed.

Risk: scatter stores dominate top-down recovery.

Mitigation: vectorize arithmetic first; leave stores scalar unless benchmark data proves scatter is the bottleneck.

## Implementation Order Summary

Recommended order:

1. Add isolated SIMD helper layer.
2. Add grouped offset caches.
3. Implement SIMD `solve_2x2_batch(...)` with scalar fallback.
4. Benchmark `2x2` dense solve impact.
5. SIMD structured `2x2` assembly for leaf, one-child, and two-child groups.
6. SIMD `recover_2x2_top_down_group(...)`.
7. Benchmark bottom-up, top-down, and total tree solve time.
8. Implement SIMD `solve_3x3_batch(...)` only if relevant cases justify it.
9. Consider packed SoA buffers if gather overhead is too high.
10. Consider SIMD counters after the path is stable.

## Success Definition

The feature is successful when:

- The serial structured tree solver uses explicit SIMD chunks for common `2x2` grouped paths.
- Scalar fallback remains available and correct.
- Existing reduced-lung tests pass.
- Large structured tree benchmarks show improved `tree_dense_s`, `tree_bottom_up_s`, `tree_top_down_s`, or total `tree_solve_s`.
- Small cases do not regress meaningfully because scalar thresholds and tails bypass SIMD overhead.
