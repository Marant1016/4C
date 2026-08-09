# Assemply Step 9.4

## Goal

Batch Kelvin-Voigt airway 3x3 dynamic tree coefficient writes without changing the assembled coefficient values.

This step touches only airway structured tree-linearization assembly. It does not change the Kelvin-Voigt equations, sparse Jacobian assembly, or tree solver logic.

## Target

File:

```text
src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.cpp
```

Function:

```text
evaluate_tree_linearization_kelvin_voigt_wall(...)
```

## What Changed

Before this step, each Kelvin-Voigt airway element performed four scalar dynamic tree coefficient replacements:

```text
momentum row, q1
momentum row, q2
mass row, q1
mass row, q2
```

After this step, the evaluator now:

```text
1. Reuses the existing resistance derivative scratch buffers for the two momentum coefficient batches.
2. Adds the inertia derivative into those buffers in one element loop.
3. Fills a persistent mass-row scratch vector with local_row_id + 1.
4. Calls replace_values(...) once for each coefficient group.
```

The four batched replacements are:

```text
target.replace_values(data.local_row_id, data.lid_q1, momentum q1 coefficients)
target.replace_values(data.local_row_id, data.lid_q2, momentum q2 coefficients)
target.replace_values(mass_row_id, data.lid_q1, mass q1 coefficients)
target.replace_values(mass_row_id, data.lid_q2, mass q2 coefficients)
```

## Why

Kelvin-Voigt airway elements contribute 3x3 structured tree blocks. The dynamic tree assembly path already had the values available in contiguous scratch buffers, but still wrote them with per-element scalar `replace_value(...)` calls.

Using `replace_values(...)` reduces per-coefficient virtual dispatch in the direct structured tree assembly path while preserving the same row, dof, and value triples.

## Expected Effect

Expected improvement is small to medium for trees with many Kelvin-Voigt airways.

The gen16 rigid-airway input is not expected to show a meaningful timing change from this step because it does not use Kelvin-Voigt airway wall mechanics heavily.

## Verification

Build the reduced-lung unit target:

```bash
cmake --build /scratch/Rodriguez/workspace/4C/4C/build/release \
  --target unittests_reduced_lung \
  --parallel 4
```

Run reduced-lung unit tests:

```bash
ctest -R "^unittests_reduced_lung$" --output-on-failure
```

## Safety

Step 9.1 added direct-vs-generic structured assembly coverage for Kelvin-Voigt airway paths. This protects the batched write change by checking that direct structured coefficients still match the generic tree linearization path and that the direct tree correction still matches the sparse solver.
