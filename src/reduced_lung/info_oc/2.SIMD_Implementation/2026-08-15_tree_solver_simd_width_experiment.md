# Tree Solver SIMD Width Experiment

Date: 2026-08-15

## Purpose

Record the temporary SIMD-width experiment for the reduced-lung tree solver.

The goal was to compare the runtime impact of different effective values of `tree_solver_simd::width()`.

## Measured Results

Input case:

```text
files/gen16_inputs/reduced_lung_lung_tree_gen16_500steps_newton_tree.4C.yaml
```

Measured `Calculation` times:

| Build | Expected width | Calculation time |
|---|---:|---:|
| `build/release-scalar` | 1 | 16.42 s |
| `build/release` | 2 | 14.86 s |
| `build/release-avx2` | 4 | 15.16 s |
| `build/release-native` | 8 | 19.24 s |

Conclusion:

- Width 2 was the fastest tested case.
- Width 1 was slower than width 2, so the default SIMD path is useful.
- Width 4 and width 8 were slower, likely because this workload is not purely arithmetic-limited.
- Wider SIMD may increase gather/scatter overhead, register pressure, masking/padding cost, and CPU frequency penalties.

## How Width Values Were Tested

### Width 1

Width 1 was forced with a temporary compile-time switch added in:

```text
src/reduced_lung/src/solver/tree/4C_reduced_lung_tree_linear_solver.cpp
```

Temporary macro:

```text
FOUR_C_REDUCED_LUNG_FORCE_SCALAR_TREE_SOLVER=1
```

Build flag:

```text
-DFOUR_C_REDUCED_LUNG_FORCE_SCALAR_TREE_SOLVER=1
```

This disables the `<experimental/simd>` path and makes `tree_solver_simd::width()` return `1`.

### Width 2

Width 2 was tested with the normal release build:

```text
build/release
```

Typical flags:

```text
-O3 -DNDEBUG
```

This uses `std::experimental::native_simd<double>` with the default compiler target, expected to be width 2 on this machine.

### Width 4

Width 4 was tested with an AVX2 build:

```text
build/release-avx2
```

Important flags:

```text
-O3 -DNDEBUG -mavx2 -mfma
```

This makes `std::experimental::native_simd<double>` use the AVX2 target, expected to hold 4 doubles per SIMD vector.

### Width 8

Width 8 was tested with a native build:

```text
build/release-native
```

Important flags:

```text
-O3 -DNDEBUG -march=native
```

On the tested CPU, `-march=native` enables AVX-512, expected to hold 8 doubles per SIMD vector.

## Important Note For Final Code

The scalar switch was added only for this experiment and should not remain in the final code unless the team wants a permanent profiling option.

Temporary code to remove later:

```text
FOUR_C_REDUCED_LUNG_FORCE_SCALAR_TREE_SOLVER
```

Reason to remove:

- It can accidentally disable SIMD.
- It was only needed to measure the width 1 case.
- The experiment already showed that width 1 is slower than width 2 for this benchmark.

## Short Summary

The experiment shows that larger SIMD width is not automatically faster for this tree-solver workload. The default release build, expected width 2, gave the best measured calculation time.
