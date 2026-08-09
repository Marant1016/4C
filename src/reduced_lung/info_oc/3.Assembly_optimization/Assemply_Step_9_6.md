# Assemply Step 9.6

## Goal

Batch Four-Element Maxwell terminal-unit dynamic tree coefficient writes without changing the assembled coefficient values.

This step touches only terminal-unit rheology tree-linearization assembly. It does not change residual assembly, sparse Jacobian assembly, elasticity evaluation, or tree solver logic.

## Target

Files:

```text
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.hpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp
```

Function:

```text
evaluate_four_element_maxwell_tree_linearization(...)
```

## What Changed

Before this step, each Four-Element Maxwell terminal unit computed one dynamic `q` coefficient and immediately wrote it with scalar `replace_value(...)`:

```text
target.replace_value(data.local_row_id[i], data.lid_q[i], grad_q)
```

After this step, the evaluator now:

```text
1. Stores dynamic q coefficients in a persistent FourElementMaxwell scratch vector.
2. Hoists references to viscosity, Maxwell elasticity, Maxwell viscosity, and reference-volume vectors.
3. Computes the dt-scaled Maxwell elasticity term once per element.
4. Calls target.replace_values(data.local_row_id, data.lid_q, grad_q) once for the full block.
```

The batched coefficient value remains:

```text
grad_q = -elastic_pressure_grad_dp_el
         - (eta + E_m * dt * eta_m / (E_m * dt + eta_m)) / v0
```

## Why

The direct structured tree assembly path can replace a contiguous group of row/dof/value triples with one `replace_values(...)` call. Four-Element Maxwell terminal units already share the same row and dof vector shape as Kelvin-Voigt terminal units, so batching removes avoidable scalar coefficient dispatch while preserving the same values.

## Expected Effect

Expected improvement is small to medium for inputs with many Four-Element Maxwell terminal units.

The gen16 rigid/linear baseline is not expected to show a meaningful timing change unless Four-Element Maxwell terminal units are active.

## Verification

Built the reduced-lung unit target:

```bash
cmake --build /scratch/Rodriguez/workspace/4C/4C/build/release \
  --target unittests_reduced_lung \
  --parallel 4
```

Ran reduced-lung unit tests:

```bash
ctest -R "^unittests_reduced_lung$" --output-on-failure
```

Result:

```text
100% tests passed, 0 tests failed out of 1
```

## Safety

Step 9.1 added direct-vs-generic structured assembly coverage for Four-Element Maxwell terminal-unit paths. That coverage protects this batched write change by checking that direct structured coefficients still match the generic tree linearization path and that the direct tree correction still matches the sparse solver.
