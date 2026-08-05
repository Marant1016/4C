# Assemply Step 5

## Goal

Separate static structured tree-linearization coefficients from dynamic coefficients.

Static row entries are now built once for a fixed `TreeLinearization` shape. Each Newton tree assembly then updates only dynamic derivative coefficients.

## Motivation

Many structured tree coefficients are constants for a fixed topology and model selection:

```text
Rigid airway pressure coefficients: +1, -1
Terminal-unit pressure coefficients: +1, -1
Junction pressure and flow coefficients: +1, -1
Boundary-condition coefficient: 1.0
```

Before this step, those constants were appended every tree-linearization assembly. That rebuilt rows even though only a subset of values changes during Newton iterations.

## Code Changes

### Dynamic Replacement API

Files:

```text
src/reduced_lung/src/4C_reduced_lung_tree_linearization.hpp
src/reduced_lung/src/4C_reduced_lung_tree_linearization.cpp
```

Added:

```cpp
void TreeLinearization::replace_value(int local_row_id, int local_dof_id, double value);
```

`replace_value(...)` updates an existing row entry and asserts that the static row pattern already contains the requested dof. It does not append new entries.

### Static Assembly Pipeline

Files:

```text
src/reduced_lung/src/4C_reduced_lung_helpers.hpp
src/reduced_lung/src/4C_reduced_lung_helpers.cpp
```

Added one-time static tree-linearization callbacks:

```cpp
using StaticTreeLinearizationAssembler =
    std::function<void(TreeLinearization& linearization)>;
```

and:

```cpp
std::vector<NamedStaticTreeLinearizationAssembler>
    tree_linearization_static_assemblers;
```

The default pipeline now registers static assemblers for:

```text
Airways
TerminalUnits
Junctions
BoundaryConditions
```

The dynamic per-iteration assembler list now contains only:

```text
Airways
TerminalUnits
```

because junction and boundary-condition tree-linearization coefficients are static.

### Newton Solver Static Initialization Guard

Files:

```text
src/reduced_lung/src/4C_reduced_lung_newton_solver.hpp
src/reduced_lung/src/4C_reduced_lung_newton_solver.cpp
```

Added:

```cpp
bool tree_linearization_static_initialized_ = false;
```

`NewtonSolver::assemble_tree_linearization_for_current_state()` now:

```text
1. Resets only when the TreeLinearization dimensions change.
2. Initializes row capacities once after reset.
3. Appends static row entries once after reset.
4. Runs dynamic assemblers every Newton tree assembly.
5. Does not clear row entries once the static pattern is initialized.
```

### Airways Static And Dynamic Split

Files:

```text
src/reduced_lung/src/airways/4C_reduced_lung_airways_common.hpp
src/reduced_lung/src/airways/4C_reduced_lung_airways.hpp
src/reduced_lung/src/airways/4C_reduced_lung_airways.cpp
src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.hpp
src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.cpp
```

Added a static tree-linearization evaluator for each airway model block.

Rigid airway static row pattern:

```text
p1: +1
p2: -1
q1: 0.0 placeholder
```

Rigid airway dynamic update:

```text
q1: -resistance_derivative - inertia_derivative
```

Kelvin-Voigt airway static row pattern:

```text
momentum row: p1 +1, p2 -1, q1 0.0, q2 0.0
mass row:     p1 +1, p2 +1, q1 0.0, q2 0.0
```

Kelvin-Voigt dynamic update replaces only q coefficients.

### Terminal-Unit Static And Dynamic Split

Files:

```text
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_common.hpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit.hpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit.cpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.hpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp
```

Terminal-unit static row pattern:

```text
p1: +1
p2: -1
q:  0.0 placeholder
```

Terminal-unit dynamic update replaces only the q derivative coefficient.

### Static-Only Junction And Boundary Contributions

Files:

```text
src/reduced_lung/src/4C_reduced_lung_junctions.cpp
src/reduced_lung/src/4C_reduced_lung_boundary_conditions.cpp
```

The existing append-only junction and boundary tree-linearization functions are now called only from the static assembly pipeline.

### Tests And Benchmarks

Updated manual `TreeLinearization` assembly helpers to run:

```text
capacity initializers
static assemblers
```

Files:

```text
src/reduced_lung/tests/4C_reduced_lung_tree_linear_solver_test.cpp
src/reduced_lung/tests/4C_reduced_lung_distributed_tree_linear_solver_test.np2.cpp
src/reduced_lung/benchmark_tests/4C_reduced_lung_solver_benchmark.cpp
src/reduced_lung/benchmark_tests/4C_reduced_lung_distributed_tree_solver_benchmark.cpp
```

## Expected Counter Impact

Expected lower repeated coefficient insertion time in:

```text
tree_assembly_airways_s
tree_assembly_terminal_units_s
tree_assembly_junctions_s
tree_assembly_boundary_conditions_s
```

`tree_assembly_junctions_s` and `tree_assembly_boundary_conditions_s` should become near zero after the first static initialization because those coefficients no longer run in the per-iteration dynamic list.

This step is not expected to directly remove:

```text
tree_assembly_solver_update_s
```

because `TreeNewtonLinearSolver::set_tree_linearization(...)` still resolves all solver coefficient values from the `TreeLinearization` rows each Newton tree assembly.

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

## Re-Measurement

From:

```text
/scratch/Rodriguez/workspace/4C/4C
```

run:

```bash
FOUR_C_REDUCED_LUNG_TREE_PROFILE=1 ./build/release/4C \
  ../files/reduced_lung_lung_tree_gen16_500steps_newton_tree.4C.yaml \
  ../output/gen16_newton_tree_profile
```

Compare against the Step 4 profile, especially the model callback counters and total tree assembly time.
