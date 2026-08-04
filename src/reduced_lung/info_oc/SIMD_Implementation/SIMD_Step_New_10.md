# SIMD Step New 10

Date: 2026-08-04

## Goal

Update reduced-lung tree linear solver tests to cover the full serial structured batch path added across the SIMD implementation steps.

The new tests exercise forced grouped traversal, strict dense-fallback policy, larger groups, odd group sizes, mixed block sizes, and solver reuse.

## Files Changed

```text
src/reduced_lung/tests/4C_reduced_lung_tree_linear_solver_test.cpp
```

## New Test Inputs

Added a generated asymmetric tree topology with 13 outlet leaves.

The topology includes:

- Leaf elements, covering `child_count == 0`.
- A serial trunk element, covering `child_count == 1`.
- Interior bifurcation elements, covering `child_count == 2`.
- An odd leaf group count larger than common SIMD widths.

Added parameter builders for:

```text
make_large_asymmetric_airway_parameters(...)
make_large_mixed_airway_terminal_unit_parameters(...)
```

The airway builder can create all-rigid `2x2` trees or all-Kelvin-Voigt `3x3` trees.

The mixed builder creates airway/terminal-unit trees with both rigid `2x2` airways, Kelvin-Voigt `3x3` airways, and terminal-unit leaves.

## New Test Helpers

Added `compare_forced_batch_structured_tree_and_sparse_corrections(...)`.

This helper:

- Compares forced-batch structured tree corrections against sparse linear solver corrections.
- Enables `force_batch_tree_solve`.
- Enables `error_on_dense_fallback`.
- Uses structured tree linearization.
- Checks forced-batch profile counters.

Added `compare_reused_forced_batch_structured_tree_solver_corrections(...)`.

This helper reuses one forced-batch structured tree solver across multiple states and compares each solve against the sparse solver.

## New Tests

Added:

```text
ForcedBatchLargeRigidAirwaysMatchSparseSolver
ForcedBatchLargeKelvinVoigtAirwaysMatchSparseSolver
ReusedForcedBatchMixedTreeSolverMatchesSparseSolver
```

Coverage summary:

- `2x2 + 0/1/2 child`: large rigid airway tree.
- `3x3 + 0/1/2 child`: large Kelvin-Voigt airway tree.
- Mixed `2x2` and `3x3` groups: mixed airway/terminal-unit tree.
- Group sizes larger than common SIMD widths: 13-leaf generated trees.
- Group sizes not divisible by common SIMD widths: 13-leaf generated trees.
- Reused structured tree solver: mixed forced-batch reuse test.

## Profile Assertions

The forced-batch tests assert:

```text
dense_solve_count == element_count * solve_count
dense_fallback_count == 0
unsupported_block_fallback_count == 0
```

When SIMD counters are nonzero, they also assert:

```text
simd_lane_count > 0
scalar_group_count == 0
scalar_tail_lane_count == 0
```

These assertions verify that normal validation inputs do not rely on dense fallback or unsupported block fallback, while keeping tests portable for non-SIMD builds.

## Behavior Preserved

- Existing small-tree tests remain unchanged.
- Existing sparse/structured correction comparisons remain unchanged.
- Existing workflow tests remain unchanged.
- The distributed tree solver tests were not changed.

## Verification

Commands run after implementation:

```text
cmake --build build/debug --target unittests_reduced_lung --parallel 4
git diff --check
ctest -R "^unittests_reduced_lung$" --output-on-failure
```
