# Optimization Step 1 - Serial Tree Solver Layout Cleanup

## Scope

This step applies the first low-risk optimization pass to the serial reduced-lung tree linear solver.

Only the serial `TreeNewtonLinearSolver` path was changed. The distributed `DistributedTreeNewtonLinearSolver` path was intentionally left unchanged because it is a separate implementation with its own communication and workspace logic.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Flattened Local Dense Matrix Storage

The serial `ElementWorkspace::matrix` storage changed from nested vectors:

```cpp
std::vector<std::vector<double>> matrix;
```

to one contiguous flat vector:

```cpp
std::vector<double> matrix;
```

The dense solver now indexes entries as:

```cpp
matrix[row * block_size + col]
```

This removes one layer of pointer indirection per matrix row and improves locality for the small dense local systems used during bottom-up condensation.

### Fixed Serial Child-Interface Storage

The serial tree solver supports at most two children per parent. The per-element child-interface storage changed from a dynamically allocated vector to fixed-size storage:

```cpp
std::array<ChildInterfacePlan, 2> child_interfaces;
int child_interface_count;
```

The matching numeric child-pressure recovery arrays in `ElementWorkspace` were also changed to fixed-size arrays:

```cpp
std::array<double, 2> child_pressure_slope;
std::array<double, 2> child_pressure_intercept;
```

This avoids per-element heap storage for child-interface metadata in the serial path.

### Split Serial Subtree Relations

The serial subtree relation storage changed from an Array-of-Structs layout:

```cpp
std::vector<SubtreeRelation> subtree_relations;
```

to separate arrays:

```cpp
std::vector<double> subtree_relation_G;
std::vector<double> subtree_relation_h;
```

This is a small Struct-of-Arrays step and prepares the bottom-up loop for later layer-wise batching or vectorization.

## Behavior Preserved

- The bottom-up/top-down tree algorithm was not changed.
- The Newton correction convention remains `J * delta = -F`.
- The structured and sparse coefficient-source behavior is unchanged.
- The serial tree solver still reuses its symbolic plan and per-element workspaces.
- The distributed tree solver implementation was not changed.

## Intended Benefit

This step reduces pointer chasing and dynamic per-element storage in the serial solver hot path. It is intended as preparation for later optimizations such as specialized `2x2`/`3x3` local solvers and direct structured coefficient access.

## Verification

Verification run successfully:

```text
git diff --check
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
```

The serial reduced-lung unit-test target passed. Benchmarks were not run in this step.
