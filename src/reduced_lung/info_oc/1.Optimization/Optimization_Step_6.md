# Optimization Step 6 - Remove Virtual Coefficient Dispatch

## Scope

This step implements Phase 6 from the serial tree solver vectorization guide.

Only the serial `TreeNewtonLinearSolver` coefficient access path was changed. The distributed `DistributedTreeNewtonLinearSolver` implementation was intentionally left unchanged.

The `3x3` batch dense solver was considered but not added in this step because it is not required for removing virtual coefficient dispatch. It remains a separate follow-up optimization.

## Files Changed

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

## Changes Made

### Removed Virtual Coefficient Provider Base

The serial coefficient access classes no longer derive from a virtual base class.

The previous path used a base pointer selected before traversal:

```cpp
const TreeCoefficientProvider* coefficients = ...;
matrix_value(*coefficients, ...);
required_matrix_value(*coefficients, ...);
```

This meant every coefficient lookup in bottom-up and top-down hot loops passed through virtual dispatch.

The serial solver now uses concrete coefficient provider objects:

```cpp
SparseTreeCoefficientProvider
StructuredTreeCoefficientProvider
```

Their `value(...)` methods are non-virtual.

### Added Compile-Time Coefficient Access

`matrix_value(...)` and `required_matrix_value(...)` are now templated on the concrete coefficient provider type:

```cpp
template <typename CoefficientProvider>
double matrix_value(const CoefficientProvider& coefficients, ...);

template <typename CoefficientProvider>
double required_matrix_value(const CoefficientProvider& coefficients, ...);
```

This lets the compiler bind coefficient access statically for the serial sparse and structured solve paths.

### Branched Once Before Traversal

The serial solve now branches once on `coefficient_source_` and invokes the same solve body with the concrete provider:

```cpp
if (coefficient_source_ == TreeNewtonLinearSolverCoefficientSource::StructuredTreeBlocks)
{
  const StructuredTreeCoefficientProvider structured_coefficients(*tree_linearization_, profile_);
  solve_with_coefficients(structured_coefficients);
}
else
{
  const SparseTreeCoefficientProvider sparse_coefficients(jacobian, profile_);
  solve_with_coefficients(sparse_coefficients);
}
```

Inside `solve_with_coefficients(...)`, bottom-up assembly, `2x2` batch dense solves, scalar dense fallbacks, and top-down recovery use the concrete provider without virtual dispatch.

## Behavior Preserved

- Sparse Jacobian and structured tree coefficient behavior are preserved.
- Coefficient lookup profiling is preserved.
- The bottom-up and top-down mathematical algorithm is unchanged.
- The Step 5 `2x2` batch solver path is unchanged.
- `3x3` blocks still use the scalar pivoted dense solver.
- No direct structured coefficient storage was added.
- `TreeLinearization::value(...)` is still used by the structured provider.
- The distributed tree solver was not changed.

## Intended Benefit

This step removes virtual dispatch from serial coefficient lookup in the hot traversal loops. It is a low-risk intermediate step before deeper structured-coefficient optimizations such as precomputed coefficient locations or direct assembly into tree-solver block storage.

## Follow-Up

The structured provider still reads from `TreeLinearization::value(...)`, which searches row entries. Removing that row-entry search is a separate follow-up and should be done with precomputed coefficient locations or direct structured block storage.

The `3x3` batch dense solver also remains a separate follow-up optimization.

## Verification

Verification passed:

```text
git diff --check
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
ctest -R "reduced_lung_.*newton_tree.*\.4C\.yaml-p1$" --output-on-failure
```

No benchmark rerun was performed in this step.
