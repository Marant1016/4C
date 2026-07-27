# Serial Tree Solver Vectorization Guide

Date: 2026-07-27

## Goal

Implement a vectorization-oriented serial tree-based solver path for `TreeNewtonLinearSolver`.

The first milestone intentionally ignores the distributed solver. Only the serial `TreeNewtonLinearSolver` should be changed until the new data layout and vectorization strategy are correct and benchmarked.

The key idea is:

```text
Tree dependencies prevent vectorizing across parent-child levels.
Elements inside the same tree layer are independent.
Therefore, vectorization should happen across elements in the same layer, grouped by compatible block shape.
```

## Current Starting Point

After `Optimization_Step_1`, the serial solver has some low-risk layout improvements:

- local dense matrices are stored as flat contiguous vectors,
- child interfaces use fixed-size arrays of length two,
- subtree relations are split into separate `G` and `h` arrays.

However, the serial solver is still mostly Array-of-Structs at the element-plan level:

```cpp
std::vector<ElementSolvePlan> element_plans_;
std::vector<ElementWorkspace> element_workspaces_;
```

This layout is easier to read but is not ideal for cache locality, batching, or compiler vectorization.

## Guiding Principles

- Prefer Structure-of-Arrays over Array-of-Structs for hot numeric data.
- Avoid irregular pointer chasing in inner loops.
- Avoid full-array writes or clears when every entry is overwritten later.
- Avoid scalar per-element loops when elements can be processed as a same-shape batch.
- Keep the mathematical algorithm unchanged until the data layout is verified.
- Keep the distributed tree solver untouched during the first milestone.

## Phase 1: Baseline Measurement

Before further changes, measure the current serial structured tree path.

Recommended benchmark entries:

```text
ReducedLung/LinearSolve/BalancedAirways/StructuredTree
ReducedLung/FullSolve/SingleTerminalUnit/NewtonTree
ReducedLung/FullSolve/SerialAirways/NewtonTree
ReducedLung/FullSolve/BalancedAirways/NewtonTree
```

Record at least:

```text
tree_solve_s
tree_bottom_up_s
tree_top_down_s
tree_dense_s
tree_lookup_s
tree_dense_solves
tree_lookups
tree_workspace_dofs
tree_max_block
```

This tells whether the dominant cost is dense local solves, coefficient lookup, top-down recovery, or general traversal overhead.

## Phase 2: Convert Serial Plan Data To SoA

Replace serial `ElementSolvePlan` hot fields with Structure-of-Arrays storage.

Target fixed per-element arrays:

```cpp
std::vector<int> global_element_id_;
std::vector<int> inlet_pressure_local_dof_;
std::vector<int> inlet_flow_unknown_index_;
std::vector<int> outlet_pressure_unknown_index_;
std::vector<int> block_size_;
std::vector<int> child_interface_count_;
std::vector<unsigned char> is_leaf_;
```

Store variable-size data in flat arrays with offsets:

```cpp
std::vector<int> unknown_offset_;
std::vector<int> equation_offset_;
std::vector<int> child_interface_offset_;

std::vector<int> unknown_global_dof_ids_;
std::vector<int> unknown_local_dof_ids_;
std::vector<int> equation_rows_;
```

Convert child-interface metadata to SoA:

```cpp
std::vector<int> child_element_index_;
std::vector<int> pressure_row_;
std::vector<int> parent_outlet_pressure_local_dof_;
std::vector<int> child_inlet_pressure_local_dof_;
std::vector<int> child_inlet_flow_local_dof_;
std::vector<int> parent_outlet_pressure_unknown_index_;
```

Keep debug strings out of the hot plan data. If detailed error messages are still needed, store debug context separately or build it only on error paths.

Definition of done for this phase:

- the scalar serial algorithm still works,
- serial reduced-lung tests pass,
- no distributed tree code is changed,
- no vectorized dense kernels are introduced yet.

## Phase 3: Build Layer And Shape Groups

Create grouped traversal data during symbolic setup.

The current traversal already has tree layers:

```cpp
tree_metadata_.bottom_up_layers
tree_metadata_.top_down_layers
```

For vectorization, build groups inside each layer by compatible shape:

```text
block_size = 2, child_count = 0
block_size = 2, child_count = 1
block_size = 2, child_count = 2
block_size = 3, child_count = 0
block_size = 3, child_count = 1
block_size = 3, child_count = 2
```

Use one flat index array plus group ranges:

```cpp
struct ElementGroup
{
  int begin = 0;
  int end = 0;
  int block_size = 0;
  int child_count = 0;
};

std::vector<int> grouped_element_indices_;
std::vector<std::vector<ElementGroup>> bottom_up_layer_groups_;
std::vector<std::vector<ElementGroup>> top_down_layer_groups_;
```

This allows the solve phase to process same-shape elements in contiguous batches.

## Phase 4: Avoid Big Array Clears

Review hot-path full-array writes such as:

```cpp
std::fill(subtree_relation_G_.begin(), subtree_relation_G_.end(), 0.0);
std::fill(subtree_relation_h_.begin(), subtree_relation_h_.end(), 0.0);
std::fill(inlet_pressure_by_element_.begin(), inlet_pressure_by_element_.end(), NaN);
delta.put_scalar(0.0);
```

If every entry is overwritten during a correct solve, avoid clearing it first.

For safety during development, use debug-only validation or stamp arrays instead of repeatedly clearing large numeric arrays:

```cpp
std::vector<int> relation_stamp_;
int current_solve_stamp_ = 0;
```

Potential rule:

```text
Release hot path: do not clear arrays that are guaranteed to be fully overwritten.
Debug path: validate that required entries were written before use.
```

For `delta`, the serial top-down recovery should write every dof exactly once. After this is verified, `delta.put_scalar(0.0)` can be removed from the serial tree solve path.

## Phase 5: Batch Dense Local Solves

Most local systems are very small:

```text
2x2: rigid airway and terminal-unit local blocks
3x3: Kelvin-Voigt airway local blocks
```

Do not start by writing explicit SIMD intrinsics. First create compiler-vectorizable batch kernels.

Add fast batch kernels:

```cpp
solve_2x2_batch(...);
solve_3x3_batch(...);
```

These should operate on contiguous SoA numeric arrays, for example for `2x2`:

```cpp
a00[i], a01[i], a10[i], a11[i]
rhs0[i], rhs1[i]
pin_rhs0[i], pin_rhs1[i]
intercept0[i], intercept1[i]
slope0[i], slope1[i]
```

Use a safe fast path plus fallback:

```text
if determinant/pivot is safely above tolerance:
  use vectorizable direct formula
else:
  use scalar pivoting fallback for that element
```

This keeps the numerical safety of the current pivoted solver while enabling vectorization for normal cases.

## Phase 6: Remove Virtual Coefficient Lookup From Hot Loops

The current serial code reads coefficients through:

```cpp
TreeCoefficientProvider::value(...)
TreeLinearization::value(...)
```

This introduces virtual dispatch and row-entry searches in the hot path.

For runtime `NewtonTree`, the coefficient source is structured. Branch once before traversal and call a structured-specific solve path instead of using virtual calls inside the element loops.

Better follow-up:

- precompute coefficient locations during symbolic setup,
- gather coefficients into SoA numeric buffers before each batch solve,
- eventually assemble structured coefficients directly into tree-solver block storage instead of into generic `TreeLinearization` rows.

## Phase 7: Vectorize Bottom-Up Batches

After SoA and grouping are available, bottom-up traversal should work by groups:

```text
for each bottom-up layer:
  for each group in layer:
    gather coefficients for group
    substitute child relations for group
    solve 2x2 or 3x3 batch
    write subtree_relation_G/H for group
```

This still uses loops, but they should be regular contiguous loops that the compiler can vectorize.

The goal is not to eliminate all loops. The goal is to eliminate irregular scalar loops when a same-shape batch can be processed through contiguous arrays.

## Phase 8: Vectorize Top-Down Batches

Top-down recovery is also layer-parallel.

For each group:

```text
known inlet_pressure[element]
recover unknown corrections
write delta values
compute child inlet pressures
```

With SoA work arrays, recovery becomes regular arithmetic:

```cpp
value0[i] = slope0[i] * inlet_pressure[i] + intercept0[i];
value1[i] = slope1[i] * inlet_pressure[i] + intercept1[i];
```

This is a natural vectorization target.

## First Milestone

The first milestone should be deliberately limited:

```text
Convert serial TreeNewtonLinearSolver plan/workspace data to SoA while keeping the existing scalar algorithm.
```

Do not include in the first milestone:

- distributed solver refactoring,
- explicit SIMD intrinsics,
- specialized `2x2` or `3x3` kernels,
- direct block assembly replacing `TreeLinearization`.

This isolates the data-layout risk from the numerical-kernel risk.

## Suggested Milestone Sequence

1. Baseline benchmark current serial structured tree path.
2. Convert serial plan/workspace data to SoA, keeping scalar behavior.
3. Build bottom-up and top-down layer groups by block size and child count.
4. Remove avoidable big-array clears after overwrite guarantees are verified.
5. Add vectorizable `2x2` batch dense solver with scalar fallback.
6. Add vectorizable `3x3` batch dense solver with scalar fallback.
7. Remove virtual structured-coefficient lookup from hot loops.
8. Benchmark small, medium, balanced-tree, and gen10 airways-only cases.

## Verification Targets

Minimum serial verification after each milestone:

```text
cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
git diff --check
```

Performance verification after vectorization milestones:

```text
./build/debug/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/LinearSolve/BalancedAirways/StructuredTree --benchmark_min_time=0.01s
./build/debug/tests/benchmarktests_reduced_lung --benchmark_filter=ReducedLung/FullSolve/.*/NewtonTree --benchmark_min_time=0.01s
```

For thesis-quality timing, repeat with a release build.

## Notes

- SoA is a means to enable contiguous batch processing; it is not automatically faster if the algorithm remains fully scalar and irregular.
- The serial solver should be optimized and benchmarked first before copying ideas into `DistributedTreeNewtonLinearSolver`.
- The safest path is to preserve the current equations, dof ordering, row ordering, and Newton correction convention throughout this work.
