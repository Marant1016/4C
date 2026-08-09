# Assemply Step 3

## Goal

Preallocate per-row storage in `TreeLinearization` so structured tree-linearization assembly does not grow each row vector dynamically while appending known coefficient patterns.

This step keeps the same coefficient values and append order. It only reserves row capacities before the first assembly for a fixed problem size.

## Motivation

After Step 2, structured tree-linearization writes use:

```cpp
TreeLinearization::append_value(...)
```

instead of the generic row-searching `set_value(...)` path. However, each row is still stored as:

```cpp
std::vector<std::pair<int, double>>
```

If the row vector does not have enough capacity, `append_value(...)` may trigger vector growth. The structured tree problem has known row sizes, so those capacities can be reserved once and reused across Newton iterations.

## Code Changes

### New Row Capacity API

Files:

```text
src/reduced_lung/src/4C_reduced_lung_tree_linearization.hpp
src/reduced_lung/src/4C_reduced_lung_tree_linearization.cpp
```

Added:

```cpp
void TreeLinearization::reserve_row_entries(int local_row_id, int entry_count);
```

Implementation:

```cpp
void TreeLinearization::reserve_row_entries(int local_row_id, int entry_count)
{
  FOUR_C_ASSERT_ALWAYS(local_row_id >= 0 && local_row_id < num_rows_, ...);
  FOUR_C_ASSERT_ALWAYS(entry_count >= 0, ...);

  rows_[static_cast<std::size_t>(local_row_id)].reserve(
      static_cast<std::size_t>(entry_count));
}
```

This reserves capacity but does not change row size or existing values.

### Capacity Initializers In The Assembly Pipeline

File:

```text
src/reduced_lung/src/4C_reduced_lung_helpers.hpp
```

Added a new callback type:

```cpp
using TreeLinearizationCapacityInitializer =
    std::function<void(TreeLinearization& linearization)>;
```

and a new pipeline vector:

```cpp
std::vector<TreeLinearizationCapacityInitializer> tree_linearization_capacity_initializers;
```

This keeps row-capacity setup tied to the same model families that already own structured tree-linearization assembly.

### Registered Expected Row Capacities

File:

```text
src/reduced_lung/src/4C_reduced_lung_helpers.cpp
```

Registered capacity initializers for:

```text
Airways
TerminalUnits
Junctions
BoundaryConditions
```

Expected capacities:

```text
Rigid airway row: 3 coefficients
Kelvin-Voigt airway momentum row: 4 coefficients
Kelvin-Voigt airway mass row: 4 coefficients
Terminal-unit row: 3 coefficients
Connection pressure row: 2 coefficients
Connection flow row: 2 coefficients
Bifurcation pressure rows: 2 coefficients each
Bifurcation flow row: 3 coefficients
Boundary-condition row: 1 coefficient
```

These capacities match the current `append_value(...)` patterns in:

```text
src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.cpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp
src/reduced_lung/src/4C_reduced_lung_junctions.cpp
src/reduced_lung/src/4C_reduced_lung_boundary_conditions.cpp
```

### One-Time Initialization In Newton Solver

Files:

```text
src/reduced_lung/src/4C_reduced_lung_newton_solver.hpp
src/reduced_lung/src/4C_reduced_lung_newton_solver.cpp
```

Added a guard:

```cpp
bool tree_linearization_capacity_initialized_ = false;
```

The capacity initializers are called in:

```cpp
NewtonSolver::assemble_tree_linearization_for_current_state()
```

only when capacity has not yet been initialized for the current row/dof dimensions:

```cpp
if (!tree_linearization_capacity_initialized_)
{
  for (const auto& initialize_capacity :
      assembly_pipeline_.tree_linearization_capacity_initializers)
  {
    initialize_capacity(tree_linearization_);
  }
  tree_linearization_capacity_initialized_ = true;
}
```

If the `TreeLinearization` dimensions ever change and `reset(...)` is called, the guard is reset to `false` so capacities are recomputed for the new shape.

## Important Nuance

`TreeLinearization::clear_values()` already preserves row vector capacity after rows have grown once. Therefore this step mainly:

- Avoids first-growth allocations during the first tree-linearization assembly.
- Makes expected row sizes explicit.
- Protects future append-only assembly from accidental repeated row growth.

For a long fixed-size run, this step is not expected to be as impactful as optimizing `tree_assembly_solver_update_s`.

## Expected Counter Impact

Potentially affected counters:

```text
tree_assembly_clear_s
tree_assembly_airways_s
tree_assembly_terminal_units_s
tree_assembly_junctions_s
tree_assembly_boundary_conditions_s
```

The effect may be small because capacity is reused after the first assembly. This step is mostly a cleanup and preparation step for later direct structured assembly work.

This step is not expected to directly reduce:

```text
tree_assembly_solver_update_s
```

because that is dominated by:

```text
TreeNewtonLinearSolver::set_tree_linearization(...)
-> TreeNewtonLinearSolver::resolve_structured_coefficient_locations(...)
```

## How To Re-Measure

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

Compare against the previous Step 1 and Step 2 profiles.

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

## Next Decision

If the next profile still shows `tree_assembly_solver_update_s` as the dominant assembly cost, the next optimization should target `TreeNewtonLinearSolver::resolve_structured_coefficient_locations(...)`, because Step 2 and Step 3 only reduce the coefficient insertion side of `TreeLinearization` assembly.
