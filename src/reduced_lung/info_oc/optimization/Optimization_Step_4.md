# Optimization Step 4 - Avoid Serial Hot-Path Clears

## Scope

This step implements Phase 4 from the serial tree solver vectorization guide.

Only the serial `TreeNewtonLinearSolver` path was changed. The distributed `DistributedTreeNewtonLinearSolver` implementation was intentionally left unchanged.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Removed Subtree Relation Clears

The serial solve path no longer clears subtree relation arrays before bottom-up traversal:

```cpp
std::fill(subtree_relation_G_.begin(), subtree_relation_G_.end(), 0.0);
std::fill(subtree_relation_h_.begin(), subtree_relation_h_.end(), 0.0);
```

These values are overwritten once per element during bottom-up traversal before any parent consumes child relations.

### Added Symbolic Group Coverage Validation

Symbolic setup now validates that grouped bottom-up and top-down traversals visit every element exactly once.

The validation also checks that each grouped element matches its group's `block_size` and `child_count`.

This supports removing solve-time clears because stale per-element relation values cannot be used when traversal coverage is complete.

### Removed Inlet-Pressure NaN Clear

The serial solve path no longer clears inlet-pressure recovery storage before top-down traversal:

```cpp
std::fill(inlet_pressure_by_element_.begin(), inlet_pressure_by_element_.end(), NaN);
```

Instead, inlet-pressure recovery now uses solve stamps:

```cpp
std::vector<int> inlet_pressure_stamp_;
int current_solve_stamp_ = 0;
```

Each written inlet pressure receives the current solve stamp. Top-down recovery asserts the stamp before reading, so stale values from previous solves are rejected without clearing the full array every solve.

### Removed Serial Delta Clear

The serial tree solve no longer calls:

```cpp
delta.put_scalar(0.0);
```

Symbolic setup now validates that element dofs cover every serial correction dof exactly once. Top-down recovery writes the inlet pressure and all element unknown corrections for every element, so the full correction vector is overwritten by the serial tree solve.

## Behavior Preserved

- The bottom-up/top-down mathematical algorithm is unchanged.
- The local dense solve is unchanged.
- The equation order and unknown order are unchanged.
- Coefficient lookup behavior is unchanged.
- Layer grouping behavior from Step 3 is unchanged.
- No vectorized kernels were added.
- No direct structured coefficient lookup path was added.
- Symbolic setup initialization remains unchanged.
- The distributed tree solver was not changed.

## Intended Benefit

This step reduces unnecessary memory traffic in the serial tree solve hot path. Arrays that are guaranteed to be overwritten are no longer cleared before the solve, and inlet-pressure correctness is checked with stamps instead of full-array NaN initialization.

This prepares the solver for later batch kernels, where avoidable memory writes would otherwise compete with grouped numeric work.

## Verification

Verification passed:

```text
git diff --check
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

No benchmark rerun was performed in this step.
