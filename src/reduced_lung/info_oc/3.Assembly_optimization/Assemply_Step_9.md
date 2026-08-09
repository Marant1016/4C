# Assemply Step 9

## Goal

Extend the Step 7&8 model-loop optimizations to the remaining airway and terminal-unit model combinations.

The direct structured assembly path is already general. This step focuses on model-specific derivative loops that were not targeted because they were not dominant in the gen16 input.

## Scope

Optimize these remaining paths:

```text
airways/4C_reduced_lung_airways_flow_resistance.cpp
- nonlinear rigid airway resistance derivative
- linear Kelvin-Voigt airway resistance derivative
- nonlinear Kelvin-Voigt airway resistance derivative
- Kelvin-Voigt inertia derivative

airways/4C_reduced_lung_airways_wall_mechanics.cpp
- Kelvin-Voigt airway 3x3 dynamic tree coefficient writes

terminal_units/4C_reduced_lung_terminal_unit_elasticity.cpp
- Ogden terminal-unit elasticity pressure gradient

terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp
- Four-element Maxwell terminal-unit tree linearization
```

## Important Rule

Do not change the equations. These optimizations should only change how values are computed and written, not the mathematical result.

## Step 1: Add Correctness Coverage First

Before optimizing a model path, add or extend direct-vs-generic structured assembly tests.

Recommended fixtures:

```text
nonlinear rigid airway
nonlinear Kelvin-Voigt airway
Ogden terminal unit
Four-element Maxwell terminal unit
mixed airway/terminal tree with these models
```

Each test should:

```text
1. Assemble the generic TreeLinearization path.
2. Assemble the direct structured coefficient path.
3. Compare every coefficient consumed by TreeNewtonLinearSolver.
4. Solve with the direct structured path.
5. Compare the correction against the sparse solver.
6. Verify tree_lookup_s / coefficient lookup count remains zero for the direct path.
```

Use the Step 8 test pattern in:

```text
src/reduced_lung/tests/4C_reduced_lung_tree_linear_solver_test.cpp
```

## Step 2: Optimize Nonlinear Rigid Airway Derivative

File:

```text
src/reduced_lung/src/airways/4C_reduced_lung_airways_flow_resistance.cpp
```

Target function:

```text
evaluate_nonlinear_flow_resistance_derivative_rigid(...)
```

Recommended changes:

```text
- Store dofs.local_values_as_span() once before the loop.
- Hoist constants using density, viscosity, and pi.
- Store ref_length[i], ref_area[i], k_turb[i], and q1 locally.
- Avoid repeated area^2 and sqrt expressions when possible.
```

Expected effect:

```text
small to medium improvement if nonlinear rigid airways are used heavily
```

## Step 3: Optimize Kelvin-Voigt Airway Resistance Derivatives

File:

```text
src/reduced_lung/src/airways/4C_reduced_lung_airways_flow_resistance.cpp
```

Target functions:

```text
evaluate_linear_flow_resistance_derivative_kelvin_voigt(...)
evaluate_nonlinear_flow_resistance_derivative_kelvin_voigt(...)
evaluate_inertia_derivative_kelvin_voigt(...)
```

Recommended changes:

```text
- Store dofs.local_values_as_span() once before each loop.
- Hoist constants such as pi, density, viscosity, and dt factors.
- Reuse local q_sum and q_difference variables.
- Reuse area, area^2, and area^3 values.
- Keep branch logic unchanged for turbulence and inertia.
```

Expected effect:

```text
medium improvement for Kelvin-Voigt airway trees, especially large 3x3 blocks
```

## Step 4: Batch Kelvin-Voigt Airway 3x3 Tree Coefficient Writes

File:

```text
src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.cpp
```

Target function:

```text
evaluate_tree_linearization_kelvin_voigt_wall(...)
```

Current behavior:

```text
Each element performs four scalar replace_value(...) calls.
```

Recommended change:

```text
1. Reuse existing scratch vectors or add persistent scratch if needed.
2. Compute batched values for:
   - momentum q1
   - momentum q2
   - mass q1
   - mass q2
3. Call replace_values(...) for each coefficient group.
```

Expected effect:

```text
reduces per-coefficient virtual dispatch for Kelvin-Voigt airway 3x3 tree assembly
```

## Step 5: Optimize Ogden Elasticity Gradient

File:

```text
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_elasticity.cpp
```

Target function:

```text
ogden_hyperelastic_pressure_gradient(...)
```

Recommended changes:

```text
- Store dofs.local_values_as_span() once before the loop.
- Hoist references to model vectors and terminal-unit data vectors.
- Store volume, reference volume, beta, kappa, and q locally.
- Avoid recomputing v0_over_vi and powers unnecessarily.
```

Expected effect:

```text
limited if std::pow dominates, but still reduces surrounding overhead
```

## Step 6: Batch Four-Element Maxwell Terminal-Unit Tree Writes

File:

```text
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp
```

Target function:

```text
evaluate_four_element_maxwell_tree_linearization(...)
```

Recommended changes:

```text
- Add persistent scratch storage to FourElementMaxwell, similar to KelvinVoigt::tree_linearization_grad_q.
- Compute all grad_q values into the scratch span.
- Use target.replace_values(data.local_row_id, data.lid_q, grad_q).
- Hoist model vector references and dt-dependent terms where possible.
```

Expected effect:

```text
small to medium improvement for Four-element Maxwell terminal-unit trees
```

## Step 7: Verify Per Model

After each model path is optimized, run:

```bash
cmake --build /scratch/Rodriguez/workspace/4C/4C/build/release \
  --target 4C benchmarktests_reduced_lung unittests_reduced_lung \
  --parallel 4
```

Then run:

```bash
ctest -R "^unittests_reduced_lung$" --output-on-failure
```

If the model has a focused benchmark fixture, run it. Otherwise, add a small benchmark before relying on timing conclusions.

## Step 8: Re-Profile Gen16 As A Regression Guard

Even though gen16 does not exercise all model paths, it should not regress.

Run:

```bash
FOUR_C_REDUCED_LUNG_TREE_PROFILE=1 ./build/release/4C \
  ../files/reduced_lung_lung_tree_gen16_500steps_newton_tree.4C.yaml \
  ../output/gen16_newton_tree_profile
```

Check:

```text
newton_total_s
tree_assembly_s
tree_assembly_airways_s
tree_assembly_terminal_units_s
tree_assembly_solver_update_s
tree_lookup_s
tree_lookups
```

Expected direct-path invariants:

```text
tree_assembly_solver_update_s: 0
tree_lookup_s: 0
tree_lookups: 0
```

## Step 9: Keep Only Permanent Improvements

Keep a change if it satisfies all of these:

```text
1. Unit tests pass.
2. Direct-vs-generic coefficient comparisons pass.
3. Sparse and direct tree corrections still match.
4. The code remains readable.
5. Timings are neutral or faster.
```

Do not keep changes that only add complexity without measurable benefit.

## Notes

The largest previous gains came from eliminating generic storage, lookup, and copy work. The remaining model-specific optimizations are expected to be smaller. Their value is mainly consistency, lower overhead for non-gen16 inputs, and better long-term protection of the structured tree fast path.
