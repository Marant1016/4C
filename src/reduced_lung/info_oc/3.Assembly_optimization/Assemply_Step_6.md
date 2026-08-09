# Assemply Step 6

## Goal

Bypass the generic `TreeLinearization` row container for the serial structured tree solver and write coefficients directly into the arrays consumed by `TreeNewtonLinearSolver`.

This step keeps the legacy `TreeLinearization` path available for distributed tree solves, sparse-compatible validation, and direct tree solver tests.

## Motivation

Before this step, serial `NewtonTree` still paid this cost every Newton linearization:

```text
model callbacks
-> TreeLinearization rows
-> TreeNewtonLinearSolver::set_tree_linearization(...)
-> TreeNewtonLinearSolver::resolve_structured_coefficient_locations(...)
-> solver coefficient value arrays
```

After Step 5, model callbacks append static coefficients once and replace dynamic coefficients only. However, the serial solver still scanned or indexed through `TreeLinearization` rows to copy values into its internal arrays.

This step changes the serial path to:

```text
model callbacks
-> TreeNewtonLinearSolver coefficient arrays
-> SIMD/batch tree solve
```

## Code Changes

### Shared Coefficient Assembly Target

File:

```text
src/reduced_lung/src/4C_reduced_lung_tree_linearization.hpp
```

Added:

```cpp
class TreeCoefficientAssemblyTarget
{
 public:
  virtual ~TreeCoefficientAssemblyTarget() = default;
  virtual void append_value(int local_row_id, int local_dof_id, double value) = 0;
  virtual void replace_value(int local_row_id, int local_dof_id, double value) = 0;
};
```

`TreeLinearization` now implements this interface, so existing row-storage assembly remains available.

### Direct Target Hook In Newton Linear Solver Interface

File:

```text
src/reduced_lung/src/4C_reduced_lung_linear_solver.hpp
```

Added:

```cpp
virtual TreeCoefficientAssemblyTarget* direct_tree_coefficient_target();
```

Default implementation returns `nullptr`.

### Serial Tree Solver Direct Target

Files:

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
```

`TreeNewtonLinearSolver` now implements `TreeCoefficientAssemblyTarget` for serial `StructuredTreeBlocks` solves.

During symbolic plan construction, it registers every solver-needed coefficient location:

```text
root_boundary_coefficient_value_
equation_inlet_pressure_coefficient_values_
matrix_coefficient_values_
child_pressure_parent_coefficient_values_
child_pressure_child_coefficient_values_
child_flow_coefficient_values_
```

The direct target stores a compact row-indexed coefficient table:

```text
local_row -> small contiguous range of {local_dof, coefficient value pointer}
```

`append_value(...)` and `replace_value(...)` now write directly to the matching slot. The row-indexed table avoids a hash lookup for every dynamic coefficient update on large trees.

### Direct Structured Provider

The serial tree solve now supports a direct provider that returns already-populated coefficient values and does not access `TreeLinearization`.

The legacy structured provider is still used when `set_tree_linearization(...)` was called, preserving old direct solver tests and fallback behavior.

### Newton Solver Direct Path

File:

```text
src/reduced_lung/src/4C_reduced_lung_newton_solver.cpp
```

`NewtonSolver::assemble_tree_linearization_for_current_state()` now checks:

```cpp
if (TreeCoefficientAssemblyTarget* direct_target =
    linear_solver_->direct_tree_coefficient_target())
```

If present, it:

```text
1. Runs static tree coefficient assemblers once against the direct target.
2. Runs dynamic tree coefficient assemblers every Newton tree assembly against the direct target.
3. Skips TreeLinearization allocation/clear/capacity setup.
4. Skips TreeNewtonLinearSolver::set_tree_linearization(...).
5. Therefore skips resolve_structured_coefficient_locations(...).
```

If no direct target is available, it falls back to the existing `TreeLinearization` path.

### Assembly Callback Generalization

Files:

```text
src/reduced_lung/src/4C_reduced_lung_helpers.hpp
src/reduced_lung/src/4C_reduced_lung_helpers.cpp
src/reduced_lung/src/airways/4C_reduced_lung_airways_common.hpp
src/reduced_lung/src/airways/4C_reduced_lung_airways.hpp
src/reduced_lung/src/airways/4C_reduced_lung_airways.cpp
src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.cpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_common.hpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit.hpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit.cpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp
src/reduced_lung/src/4C_reduced_lung_junctions.hpp
src/reduced_lung/src/4C_reduced_lung_junctions.cpp
src/reduced_lung/src/4C_reduced_lung_boundary_conditions.hpp
src/reduced_lung/src/4C_reduced_lung_boundary_conditions.cpp
```

Static and dynamic structured-tree assembly callbacks now target `TreeCoefficientAssemblyTarget&` where they only append or replace coefficients.

This lets the same model coefficient code write to either:

```text
TreeLinearization
TreeNewtonLinearSolver direct coefficient target
```

## Expected Counter Impact

The main expected improvement is removal of the serial solver-update phase:

```text
tree_assembly_solver_update_s -> near zero
tree_lookup_s -> zero for direct serial path
tree_lookups -> zero for direct serial path
```

The model assembly counters still measure dynamic coefficient evaluation and writes:

```text
tree_assembly_airways_s
tree_assembly_terminal_units_s
```

Junction and boundary-condition dynamic counters should remain zero after Step 5 because those are static-only.

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

Ran a filtered serial NewtonTree benchmark to exercise the direct path:

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

Ran the real gen16 serial NewtonTree profile:

```bash
FOUR_C_REDUCED_LUNG_TREE_PROFILE=1 ./build/release/4C \
  ../files/reduced_lung_lung_tree_gen16_500steps_newton_tree.4C.yaml \
  ../output/gen16_newton_tree_profile
```

Observed:

```text
newton_total_s: 17.0854
tree_assembly_s: 1.39912
tree_assembly_solver_update_s: 0
tree_lookup_s: 0
tree_lookups: 0
tree_solve_s: 4.70658
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

Compare against Step 5, especially:

```text
tree_assembly_solver_update_s
tree_assembly_s
newton_total_s
```
