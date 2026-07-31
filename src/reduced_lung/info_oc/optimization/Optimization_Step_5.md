# Optimization Step 5 - Serial 2x2 Batch Dense Solver

## Scope

This step implements Phase 5 from the serial tree solver vectorization guide for `2x2` local blocks only.

Only the serial `TreeNewtonLinearSolver` path was changed. The distributed `DistributedTreeNewtonLinearSolver` implementation was intentionally left unchanged.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Added 2x2 Batch Workspace

The serial solver now owns reusable SoA buffers for grouped `2x2` dense solves:

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
std::vector<int> batch_2x2_fallback_lanes_;
```

These buffers are sized once during symbolic setup from the largest bottom-up group with `block_size == 2`.

### Added Safe 2x2 Batch Solve Helper

A new `solve_2x2_batch(...)` helper gathers each grouped element's local matrix and RHS values into SoA arrays and solves safe lanes with direct `2x2` formulas:

```text
det = a00 * a11 - a01 * a10
x0 = (rhs0 * a11 - a01 * rhs1) / det
x1 = (a00 * rhs1 - rhs0 * a10) / det
```

The helper solves both RHS systems used by the tree solver:

- constant RHS, producing `intercept`
- inlet-pressure RHS, producing `slope`

If the determinant is not finite or its absolute value is not safely above `pivot_tolerance_`, that lane is solved with the existing pivoted scalar `solve_dense_system(...)` fallback.

### Used Batch Path For 2x2 Groups

Bottom-up traversal now uses the Step 3 group data:

- all elements in a group are assembled first,
- groups with `block_size == 2` call `solve_2x2_batch(...)`,
- all other groups keep the existing scalar dense solve,
- subtree relations are written after each group is solved.

The profiling counter still counts solved elements, not batch calls.

## Behavior Preserved

- The bottom-up and top-down mathematical algorithm is unchanged.
- Coefficient lookup and local assembly are still scalar.
- `3x3` and larger blocks still use the existing pivoted scalar solver.
- Unsafe `2x2` lanes still use the existing pivoted scalar solver.
- The equation order and unknown order are unchanged.
- No SIMD intrinsics were added.
- No direct structured coefficient lookup path was added.
- The distributed tree solver was not changed.

## Intended Benefit

This step creates the first compiler-vectorizable dense local solve path for same-shape layer groups. It targets the common `2x2` rigid-airway and terminal-unit local blocks while preserving the numerical safety of the existing pivoted solver through per-lane fallback.

## Verification

Verification passed:

```text
git diff --check
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

No benchmark rerun was performed in this step.
