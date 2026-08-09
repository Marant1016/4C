# Assemply Step 7&8

## Goal

Finish the structured tree assembly optimization by:

1. Batching hot dynamic model coefficient writes after Steps 1-6 removed the generic storage and coefficient-resolution bottlenecks.
2. Adding direct-vs-generic structured assembly correctness tests so future fast-path changes are guarded.

## Step 7: Batched Dynamic Model Coefficient Loops

The gen16 profile after Step 6 showed arithmetic/model assembly was still measurable:

```text
tree_assembly_s: 1.39912
tree_assembly_airways_s: 0.946144
tree_assembly_terminal_units_s: 0.447769
tree_assembly_solver_update_s: 0
```

The remaining work was no longer coefficient resolution. It was model derivative evaluation plus scalar coefficient writes.

### Bulk Coefficient Replacement

File:

```text
src/reduced_lung/src/4C_reduced_lung_tree_linearization.hpp
src/reduced_lung/src/4C_reduced_lung_tree_linearization.cpp
```

Added a batched target API:

```cpp
virtual void replace_values(std::span<const int> local_row_ids,
    std::span<const int> local_dof_ids, std::span<const double> values);
```

The default implementation loops over scalar `replace_value(...)`, preserving compatibility for all existing targets.

`TreeLinearization` overrides this API directly.

File:

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
```

`TreeNewtonLinearSolver` overrides `replace_values(...)` for the direct structured coefficient target. This removes one virtual dispatch per dynamic coefficient write and keeps direct writes in the row-indexed coefficient table from Step 6.

### Rigid Airway Linear Resistance Path

Files:

```text
src/reduced_lung/src/airways/4C_reduced_lung_airways_flow_resistance.cpp
src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.cpp
```

Changes:

```text
ComputePoiseuilleResistance:
- hoists 8*pi*mu out of the element loop
- avoids repeated area indexing for area^2

evaluate_tree_linearization_rigid_wall:
- reuses the resistance-derivative scratch buffer for final q coefficient values
- performs one batched replace_values(...) call for all rigid airway q coefficients
```

This is the most relevant airway path for the gen16 YAML.

### Terminal-Unit Linear Elasticity And Kelvin-Voigt Path

Files:

```text
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_elasticity.cpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.hpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp
```

Changes:

```text
linear_elastic_pressure_gradient:
- hoists vector references and writes directly through the gradient vector reference

KelvinVoigt rheology:
- adds persistent tree_linearization_grad_q scratch storage
- computes all terminal-unit q coefficients into the scratch span
- performs one batched replace_values(...) call for the model block
```

This is the most relevant terminal-unit path for the gen16 YAML.

## Step 8: Structured Assembly Correctness Tests

Files:

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
src/reduced_lung/tests/4C_reduced_lung_tree_linear_solver_test.cpp
```

Added a solver-consumed coefficient snapshot:

```cpp
struct TreeStructuredCoefficientValue
{
  int local_row = -1;
  int local_dof = -1;
  double value = 0.0;
  const char* context = nullptr;
};
```

and:

```cpp
std::vector<TreeStructuredCoefficientValue> structured_coefficient_values() const;
```

The test suite now checks:

```text
1. Assemble generic TreeLinearization.
2. Resolve generic values into a structured TreeNewtonLinearSolver.
3. Assemble direct structured coefficients through direct_tree_coefficient_target().
4. Compare every coefficient value consumed by TreeNewtonLinearSolver.
5. Solve the direct structured path and compare the correction against sparse.
6. Verify direct coefficient lookups remain zero.
```

New tests:

```text
DirectStructuredAssemblyRigidAirwaysMatchesGenericPath
DirectStructuredAssemblyKelvinVoigtAirwaysMatchesGenericPath
DirectStructuredAssemblyMixedTerminalUnitsMatchesGenericPath
```

The fixtures cover:

```text
rigid 2x2 airway blocks
Kelvin-Voigt 3x3 airway blocks
mixed terminal-unit blocks, including Kelvin-Voigt terminal units
```

## Verification

Built release targets:

```bash
cmake --build /scratch/Rodriguez/workspace/4C/4C/build/release \
  --target 4C benchmarktests_reduced_lung unittests_reduced_lung \
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

Ran focused serial NewtonTree benchmark:

```bash
./tests/benchmarktests_reduced_lung \
  --benchmark_filter="^ReducedLung/FullSolve/SerialAirways/NewtonTree$" \
  --benchmark_min_time=0.01 \
  --benchmark_repetitions=1
```

Observed direct-path counters:

```text
tree_assembly_solver_update_s=0
tree_lookup_s=0
tree_lookups=0
```

Ran real gen16 serial NewtonTree profile:

```bash
FOUR_C_REDUCED_LUNG_TREE_PROFILE=1 ./build/release/4C \
  ../files/reduced_lung_lung_tree_gen16_500steps_newton_tree.4C.yaml \
  ../output/gen16_newton_tree_profile
```

Observed after Step 7&8:

```text
newton_total_s: 15.9412
tree_assembly_s: 1.26883
tree_assembly_airways_s: 0.856233
tree_assembly_terminal_units_s: 0.407285
tree_assembly_solver_update_s: 0
tree_lookup_s: 0
tree_lookups: 0
tree_solve_s: 4.29166
```

Compared to Step 6:

```text
newton_total_s: 17.0854 -> 15.9412
tree_assembly_s: 1.39912 -> 1.26883
tree_assembly_airways_s: 0.946144 -> 0.856233
tree_assembly_terminal_units_s: 0.447769 -> 0.407285
tree_assembly_solver_update_s: 0 -> 0
```

## Notes

No explicit experimental SIMD was added in this step. The permanent optimization is the batched coefficient path: it makes the arithmetic loops easier for the compiler to optimize and removes per-coefficient virtual dispatch on direct structured assembly.
