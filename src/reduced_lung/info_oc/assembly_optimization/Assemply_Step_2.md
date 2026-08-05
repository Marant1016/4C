# Assemply Step 2

## Goal

Avoid row search overhead in structured tree linearization coefficient insertion. This step adds a fast append-only insertion path to `TreeLinearization` and migrates the current structured tree-linearization assemblers to use it.

This step intentionally does not change the solver algorithm, coefficient values, row ordering assumptions, or the existing generic `set_value(...)` behavior.

## Motivation

Before this step, every tree-linearization coefficient write used:

```cpp
TreeLinearization::set_value(int local_row_id, int local_dof_id, double value)
```

`set_value(...)` searches the target row with `std::find_if(...)` to support insert-or-replace semantics:

```cpp
const auto entry = std::find_if(row.begin(), row.end(),
    [local_dof_id](const auto& coefficient) { return coefficient.first == local_dof_id; });
```

For structured tree assembly, the normal pattern is:

```text
TreeLinearization::clear_values()
-> write each row's known coefficients once
-> TreeNewtonLinearSolver::set_tree_linearization(...)
```

Because the rows are cleared before assembly and each row/dof pair is written once by its owning model family, replacement search is unnecessary in the hot path.

## Code Changes

### New Append API

Files:

```text
src/reduced_lung/src/4C_reduced_lung_tree_linearization.hpp
src/reduced_lung/src/4C_reduced_lung_tree_linearization.cpp
```

Added:

```cpp
void TreeLinearization::append_value(int local_row_id, int local_dof_id, double value);
```

Implementation:

```cpp
void TreeLinearization::append_value(int local_row_id, int local_dof_id, double value)
{
  FOUR_C_ASSERT_ALWAYS(local_row_id >= 0 && local_row_id < num_rows_, ...);
  FOUR_C_ASSERT_ALWAYS(local_dof_id >= 0 && local_dof_id < num_dofs_, ...);

  rows_[static_cast<std::size_t>(local_row_id)].emplace_back(local_dof_id, value);
}
```

This keeps bounds checking but skips the per-row `std::find_if(...)` duplicate/replacement search.

### Existing Generic API Preserved

`TreeLinearization::set_value(...)` was left unchanged.

It still supports safe insert-or-replace behavior and remains available for any future generic path that needs replacement semantics.

The intended API split is now:

```text
set_value(...)    = generic insert-or-replace path
append_value(...) = fast append-only structured assembly path
```

### Migrated Structured Tree Assembly Writes

Changed structured tree-linearization coefficient writes from `set_value(...)` to `append_value(...)` in:

```text
src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.cpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp
src/reduced_lung/src/4C_reduced_lung_junctions.cpp
src/reduced_lung/src/4C_reduced_lung_boundary_conditions.cpp
```

Examples:

```cpp
target.set_value(data.local_row_id[i], data.lid_p1[i], 1.0);
```

became:

```cpp
target.append_value(data.local_row_id[i], data.lid_p1[i], 1.0);
```

and:

```cpp
linearization.set_value(flow_row, local_dof_ids[BifurcationData::q_out_parent], 1.0);
```

became:

```cpp
linearization.append_value(flow_row, local_dof_ids[BifurcationData::q_out_parent], 1.0);
```

## Why This Is Safe For The Current Structured Path

The current tree-linearization assembly flow clears rows before writing:

```text
NewtonSolver::assemble_tree_linearization_for_current_state()
-> TreeLinearization::clear_values()
-> Airways::update_tree_linearization(...)
-> TerminalUnits::update_tree_linearization(...)
-> Junctions::update_tree_linearization(...)
-> BoundaryConditions::update_tree_linearization(...)
```

Each coefficient row belongs to one equation family:

- Airway state equations.
- Terminal-unit equations.
- Junction equations.
- Boundary-condition equations.

Within each family, each row/dof coefficient is written once in a fixed pattern. Therefore append-only insertion preserves the same row entries without doing redundant replacement searches.

## Expected Counter Impact

This step should affect the model callback portions of `tree_assembly_s`:

```text
tree_assembly_airways_s
tree_assembly_terminal_units_s
tree_assembly_junctions_s
tree_assembly_boundary_conditions_s
```

It is not expected to directly reduce the dominant `tree_assembly_solver_update_s`, because that counter is mostly:

```text
TreeNewtonLinearSolver::set_tree_linearization(...)
-> TreeNewtonLinearSolver::resolve_structured_coefficient_locations(...)
```

The latest profile before this step showed `tree_assembly_solver_update_s` was the largest assembly sub-cost, so Step 2 is useful but not expected to solve the main bottleneck alone.

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

Compare these values against the previous run:

```text
tree_assembly_airways_s
tree_assembly_terminal_units_s
tree_assembly_junctions_s
tree_assembly_boundary_conditions_s
tree_assembly_s
```

The most important remaining value to watch is still:

```text
tree_assembly_solver_update_s
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

## Next Decision

After re-measuring, if the model callback counters decrease but `tree_assembly_solver_update_s` remains dominant, the next optimization should target `TreeNewtonLinearSolver::resolve_structured_coefficient_locations(...)` or a direct structured coefficient update path.
