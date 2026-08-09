# Optimization Step 7 - Vectorize Serial Bottom-Up Batches

## Scope

This step implements Phase 7 from the serial tree solver vectorization guide.

Only the serial `TreeNewtonLinearSolver` bottom-up path was changed. The distributed `DistributedTreeNewtonLinearSolver` implementation was intentionally left unchanged.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Added 3x3 Batch Dense Solver

The serial solver now has reusable SoA workspace for grouped `3x3` dense solves:

```cpp
std::vector<double> batch_3x3_a00_;
std::vector<double> batch_3x3_a01_;
std::vector<double> batch_3x3_a02_;
std::vector<double> batch_3x3_a10_;
std::vector<double> batch_3x3_a11_;
std::vector<double> batch_3x3_a12_;
std::vector<double> batch_3x3_a20_;
std::vector<double> batch_3x3_a21_;
std::vector<double> batch_3x3_a22_;
```

Additional buffers store the two RHS systems, the `intercept` and `slope` outputs, and fallback lane indices.

The new `solve_3x3_batch(...)` helper gathers each grouped element's local matrix and RHS values into SoA buffers and uses a no-pivot LU fast path with pivot checks. If any lane has an unsafe or non-finite pivot, it falls back to the existing pivoted scalar `solve_dense_system(...)` for that lane.

### Sized 3x3 Batch Workspace During Symbolic Setup

Symbolic setup now computes the largest bottom-up group with `block_size == 3` and sizes the reusable `3x3` batch buffers once.

The existing `2x2` batch workspace sizing is unchanged.

### Refactored Bottom-Up Assembly To Group Phases

Bottom-up traversal now works in explicit group phases:

```text
assemble_group(group)
solve group with 2x2, 3x3, or scalar fallback
write_subtree_relation_group(group)
```

The new `assemble_group(...)` loop processes same-shape groups with equation index as the outer loop and grouped element lanes as the inner loop. Downstream child-relation substitution is also performed as a group-level loop over the group's child slots.

### Added 3x3 Batch Dispatch

Bottom-up dense solve dispatch is now:

```text
block_size == 2: solve_2x2_batch(...)
block_size == 3: solve_3x3_batch(...)
otherwise: scalar solve_dense_system(...)
```

Dense-solve profiling still counts solved elements, not batch calls.

## Behavior Preserved

- The bottom-up and top-down mathematical algorithm is unchanged.
- Coefficient access behavior is unchanged from Step 6.
- The existing `2x2` batch solver is preserved.
- Unsafe `2x2` and `3x3` lanes still use the pivoted scalar fallback.
- Non-`2x2`/`3x3` groups still use the scalar pivoted solver.
- The equation order and unknown order are unchanged.
- No SIMD intrinsics were added.
- No top-down vectorization was added.
- No direct structured coefficient storage was added.
- The distributed tree solver was not changed.

## Intended Benefit

This step makes bottom-up traversal match the intended grouped execution model more closely. Instead of treating grouping only as an ordering change, the bottom-up phase now assembles, solves, and writes subtree relations at group granularity.

This creates regular same-shape loops that are more suitable for compiler vectorization and prepares the path for future structured-coefficient gather improvements.

## Follow-Up

Structured coefficient access still calls `TreeLinearization::value(...)`, which searches entries in each row. A future optimization should precompute coefficient locations or assemble directly into tree-solver block storage.

Top-down recovery is still scalar inside each group and remains the target of Phase 8.

## Verification

Verification passed:

```text
git diff --check
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

No benchmark rerun was performed in this step.
