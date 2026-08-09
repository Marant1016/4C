# Structured Tree Linearization Assembly Optimization Guide

Date: 2026-08-04

## Motivation

The serial structured `NewtonTree` SIMD/batch workflow is already active for the gen16 rigid-tree YAML case. The profiling output shows that the remaining dominant measured cost is now structured tree linearization assembly, not the SIMD dense solve kernel.

Observed gen16 profile:

```text
newton_solves: 500
tree_linear_solves: 1000
newton_total_s: 28.0862
tree_assembly_s: 15.9798
linear_solve_s: 4.173
tree_solve_s: 4.17204
tree_dense_s: 0.687825
tree_simd_lanes: 92459000
tree_dense_fallbacks: 0
tree_unsupported_fallbacks: 0
tree_max_block: 2
```

Interpretation:

- SIMD/batch tree solve is working: all `92,459,000` element solves were counted as SIMD/batch lanes.
- Dense solve fallback is not a problem: `tree_dense_fallbacks = 0`.
- Unsupported block fallback is not a problem: `tree_unsupported_fallbacks = 0`.
- The dense kernel is small relative to total Newton time: `tree_dense_s = 0.687825` out of `28.0862` seconds.
- The next meaningful target is `tree_assembly_s = 15.9798` seconds.

## What `tree_assembly_s` Measures

`tree_assembly_s` is accumulated in:

```text
src/reduced_lung/src/4C_reduced_lung_newton_solver.cpp
NewtonSolver::assemble_tree_linearization_for_current_state()
```

It includes:

1. Clearing or resetting `TreeLinearization`.
2. Running all `assembly_pipeline_.tree_linearization_assemblers`.
3. Calling `linear_solver_->set_tree_linearization(tree_linearization_)`.
4. For the serial structured tree solver, resolving cached structured coefficient values.

It does not include:

- `TreeNewtonLinearSolver::solve(...)` bottom-up traversal.
- `solve_2x2_batch(...)` or `solve_3x3_batch(...)`.
- Top-down recovery.
- Dense fallback time.

## Main Code Paths

### Assembly Driver

```text
src/reduced_lung/src/4C_reduced_lung_newton_solver.cpp
NewtonSolver::assemble_tree_linearization_for_current_state()
```

Important operations:

```cpp
tree_linearization_.clear_values();

for (const auto& assemble_tree_linearization : assembly_pipeline_.tree_linearization_assemblers)
{
  assemble_tree_linearization(tree_linearization_, locally_relevant_dofs_, current_time_, dt_);
}

linear_solver_->set_tree_linearization(tree_linearization_);
```

### Pipeline Registration

```text
src/reduced_lung/src/4C_reduced_lung_helpers.cpp
create_default_reduced_lung_assembly_pipeline(...)
```

Registered structured tree assemblers:

```cpp
Airways::update_tree_linearization(...);
TerminalUnits::update_tree_linearization(...);
Junctions::update_tree_linearization(...);
BoundaryConditions::update_tree_linearization(...);
```

### TreeLinearization Storage

```text
src/reduced_lung/src/4C_reduced_lung_tree_linearization.cpp
```

Current storage API:

```cpp
TreeLinearization::clear_values();
TreeLinearization::set_value(int local_row_id, int local_dof_id, double value);
TreeLinearization::entries(int local_row_id) const;
```

`set_value(...)` currently searches each row with `std::find_if` before inserting or replacing a coefficient.

### Airway Tree Linearization

```text
src/reduced_lung/src/airways/4C_reduced_lung_airways.cpp
Airways::update_tree_linearization(...)
```

For the gen16 rigid airway case:

```text
src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.cpp
evaluate_tree_linearization_rigid_wall(...)
```

This loops over airway elements and writes three coefficients per rigid airway row:

```cpp
target.set_value(data.local_row_id[i], data.lid_p1[i], 1.0);
target.set_value(data.local_row_id[i], data.lid_p2[i], -1.0);
target.set_value(data.local_row_id[i], data.lid_q1[i], -resistance_derivative[i] - inertia_derivative[i]);
```

### Terminal-Unit Tree Linearization

```text
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit.cpp
TerminalUnits::update_tree_linearization(...)
```

For the gen16 terminal-unit case:

```text
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_elasticity.cpp
linear_elastic_pressure_gradient(...)

src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp
evaluate_kelvin_voigt_tree_linearization(...)
```

This loops over terminal units and writes three coefficients per terminal unit row:

```cpp
target.set_value(data.local_row_id[i], data.lid_p1[i], 1.0);
target.set_value(data.local_row_id[i], data.lid_p2[i], -1.0);
target.set_value(data.local_row_id[i], data.lid_q[i], grad_q);
```

### Junction And Boundary Linearization

```text
src/reduced_lung/src/4C_reduced_lung_junctions.cpp
Junctions::update_tree_linearization(...)

src/reduced_lung/src/4C_reduced_lung_boundary_conditions.cpp
BoundaryConditions::update_tree_linearization(...)
```

These write pressure-continuity, flow-conservation, and boundary-condition coefficients through the same `TreeLinearization::set_value(...)` path.

### Structured Coefficient Resolve

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
TreeNewtonLinearSolver::set_tree_linearization(...)
TreeNewtonLinearSolver::resolve_structured_coefficient_locations(...)
```

This copies coefficient values from `TreeLinearization` rows into solver arrays such as:

```cpp
equation_inlet_pressure_coefficient_values_
matrix_coefficient_values_
child_pressure_parent_coefficient_values_
child_pressure_child_coefficient_values_
child_flow_coefficient_values_
```

## Optimization Goal

Reduce `tree_assembly_s` for serial structured `NewtonTree`, especially for large rigid `2x2` trees, without weakening the sparse validation path or the existing structured solver correctness.

The key idea is to avoid rebuilding and searching a generic row-wise coefficient container when the serial tree solver already knows the exact symbolic coefficient locations and the coefficient shapes are fixed.

## Recommended Step 1: Add More Fine-Grained Timing

Before changing data structures, split `tree_assembly_s` into sub-counters.

Add counters to `NewtonSolverProfile` or a new structured assembly profile:

```cpp
tree_linearization_clear_time
tree_linearization_airway_time
tree_linearization_terminal_time
tree_linearization_junction_time
tree_linearization_boundary_time
tree_linearization_resolve_time
```

Minimal implementation options:

- Add explicit timing around each tree-linearization assembler callback in `NewtonSolver::assemble_tree_linearization_for_current_state()`.
- Give each callback a name in the pipeline, or register tree assemblers as a small struct `{name, callback}`.
- Alternatively add scoped timers inside each update function for a temporary investigation branch.

Expected value:

- This tells whether the main cost is coefficient calculation, `TreeLinearization::set_value(...)`, `clear_values()`, or `resolve_structured_coefficient_locations(...)`.

## Recommended Step 2: Avoid Row Search In `TreeLinearization::set_value(...)`

Current issue:

```cpp
const auto entry = std::find_if(row.begin(), row.end(), ...);
```

Every coefficient insertion checks whether the coefficient already exists. For the structured tree assembly, most rows have a small fixed pattern and are rebuilt every Newton iteration. This makes the generic search unnecessary in the common path.

Possible improvements:

1. Add `append_value(...)` for assembly after `clear_values()`.
2. Use `set_value(...)` only for paths that genuinely need replacement semantics.
3. Change tree-linearization assemblers to call `append_value(...)` when each row is known to be written once.

Example API:

```cpp
void TreeLinearization::append_value(int local_row_id, int local_dof_id, double value)
{
  rows_[static_cast<std::size_t>(local_row_id)].emplace_back(local_dof_id, value);
}
```

Validation:

- In assertion builds, optionally check duplicates after assembly or inside `append_value(...)` only when expensive checks are enabled.
- Keep `set_value(...)` unchanged for safety where needed.

Expected impact:

- Lower coefficient insertion overhead across airways, terminal units, junctions, and boundary conditions.
- This is a small, low-risk first optimization.

## Recommended Step 3: Preallocate Row Capacities

Current issue:

`TreeLinearization::clear_values()` clears every row, but rows may still grow dynamically if capacity is insufficient.

For the structured tree problem, row nonzero counts are known from metadata:

- Rigid airway row: usually 3 coefficients.
- Terminal-unit row: usually 3 coefficients.
- Connection pressure row: 2 coefficients.
- Connection flow row: 2 coefficients.
- Bifurcation pressure rows: 2 coefficients each.
- Bifurcation flow row: 3 coefficients.
- Boundary row: 1 coefficient.

Possible improvements:

1. Add a one-time `reserve_row_entries(row, count)` or `reset_with_row_capacities(...)` API.
2. During setup, compute expected row capacities from models or tree metadata.
3. Reserve once and reuse across Newton iterations.

Expected impact:

- Avoid repeated vector allocations or growth during assembly.
- Helps large cases with many rows, such as gen16.

## Recommended Step 4: Avoid Temporary Allocations In Model Coefficient Evaluators

Current issue:

Some tree-linearization callbacks create temporary `std::vector<double>` objects for derivative values, then immediately consume those vectors to write coefficients into `TreeLinearization`. For large trees and many Newton linear solves, this can add repeated heap allocation, copy/move, and cache pressure.

Important candidate paths:

```text
src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.cpp
WallMechanics::make_tree_linearization_evaluator(...)

src/reduced_lung/src/airways/4C_reduced_lung_airways_flow_resistance.cpp
FlowResistance derivative evaluators
```

The rigid airway tree-linearization path currently creates derivative vectors before writing the final row coefficients:

```cpp
auto resistance_derivative =
    resistance_derivative_evaluator(airway_data, locally_relevant_dofs, dt);
auto inertia_derivative = ...;
evaluate_tree_linearization_rigid_wall(
    target, airway_data, resistance_derivative, inertia_derivative);
```

The gen16 rigid case is especially worth checking because it has `61,639` airways and `1000` tree linear solves in the profiled run.

Terminal-unit linear elasticity is already better in this respect because `linear_elastic_pressure_gradient(...)` returns a reference to model-owned storage:

```cpp
std::vector<double>& linear_elastic_pressure_gradient(...)
```

Still, terminal-unit and Kelvin-Voigt paths should be reviewed for avoidable temporaries before adding SIMD.

Possible improvements:

1. Precompute constant rigid-linear airway derivatives once when the model is created, if they do not depend on the current Newton state.
2. Store reusable scratch arrays in airway model data for derivatives that must be recomputed each iteration.
3. Change evaluator APIs from `return std::vector<double>` to `write_into(std::span<double>)` or `write_into(std::vector<double>&)`.
4. Fuse derivative computation and coefficient insertion for simple cases, so no intermediate derivative vector is needed.
5. In the future direct structured assembly path, write dynamic coefficients directly into solver coefficient arrays and skip both temporary vectors and `TreeLinearization::set_value(...)`.

Expected impact:

- Reduces heap allocation during every tree linearization.
- Improves cache locality in the hot assembly path.
- Makes later SIMD/batch coefficient evaluation cleaner because output buffers are explicit and reusable.

## Recommended Step 5: Separate Static And Dynamic Coefficients

Many structured tree coefficients are constants for a fixed topology and model type.

Examples:

- Rigid airway pressure coefficients: `+1`, `-1`.
- Terminal-unit pressure coefficients: `+1`, `-1`.
- Junction pressure-continuity and flow-conservation coefficients: mostly `+1`, `-1`.
- Boundary-condition coefficient: `1.0`.

Only a subset is dynamic:

- Rigid airway flow resistance/inertia derivative coefficient.
- Terminal-unit `grad_q` coefficient.
- Kelvin-Voigt or nonlinear model derivative coefficients when enabled.

Possible improvements:

1. Build the static part of `TreeLinearization` once during setup.
2. During each Newton iteration, update only dynamic coefficient values.
3. Keep row entry order stable so `TreeNewtonLinearSolver::resolve_structured_coefficient_locations(...)` can keep valid `structured_entry_index` values.

Expected impact:

- Reduces repeated writes of constant coefficients.
- Reduces duplicate row construction work.
- Makes later direct structured coefficient updates easier.

## Recommended Step 6: Direct Structured Coefficient Assembly For Serial Tree Solver

Longer-term target:

Bypass the generic `TreeLinearization` row container for the serial structured tree solver. Instead, write dynamic coefficients directly into the arrays already consumed by `TreeNewtonLinearSolver`.

Current flow:

```text
model callbacks
-> TreeLinearization::set_value(...)
-> TreeNewtonLinearSolver::set_tree_linearization(...)
-> resolve_structured_coefficient_locations(...)
-> solver coefficient value arrays
-> SIMD/batch tree solve
```

Proposed fast flow:

```text
model callbacks or dedicated structured assembler
-> solver coefficient value arrays
-> SIMD/batch tree solve
```

Possible design:

1. Add a `StructuredTreeCoefficientValues` object owned by the Newton solver or tree solver.
2. Build symbolic locations once from tree metadata.
3. Add dedicated update callbacks for dynamic airway and terminal-unit coefficients.
4. Keep the old `TreeLinearization` path for validation and for distributed/sparse-compatible paths.

Expected impact:

- Removes generic row insertion cost.
- Removes `resolve_structured_coefficient_locations(...)` from every Newton linearization.
- Aligns assembly with the SoA layout used by the SIMD/batch tree solver.

Risk:

- Higher implementation complexity.
- Must preserve sparse-vs-structured validation tests.
- Must ensure distributed tree solver requirements are not broken.

## Recommended Step 7: Batch/SIMD The Dynamic Model Coefficient Loops

After the data path is cleaned up, consider SIMD for the model derivative loops themselves.

Candidate loops:

```text
airways/4C_reduced_lung_airways_flow_resistance.cpp
ComputePoiseuilleResistance and derivative evaluators

terminal_units/4C_reduced_lung_terminal_unit_elasticity.cpp
linear_elastic_pressure_gradient(...)

terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp
evaluate_kelvin_voigt_tree_linearization(...)
```

For the gen16 YAML, the most relevant cases are:

- Rigid airway linear resistance derivative.
- Terminal-unit linear elasticity pressure gradient.
- Terminal-unit Kelvin-Voigt tree linearization.

Do this only after Steps 1-6 show arithmetic is still significant. The current code likely spends substantial time in generic storage, temporary allocation, and coefficient resolve rather than floating-point arithmetic.

## Recommended Step 8: Add A Structured Assembly Correctness Test

Any fast assembly path should be tested against the existing generic path.

Suggested tests:

1. Assemble `TreeLinearization` with the existing generic path.
2. Assemble direct structured coefficient arrays with the new path.
3. Compare all coefficient values consumed by `TreeNewtonLinearSolver`.
4. Reuse rigid `2x2` and Kelvin-Voigt `3x3` test fixtures.
5. Verify the final correction still matches the sparse solver.

Existing relevant tests are in:

```text
src/reduced_lung/tests/4C_reduced_lung_tree_linear_solver_test.cpp
```

## Proposed Implementation Order

1. Add fine-grained profile counters for tree linearization assembly phases.
2. Add `TreeLinearization::append_value(...)` and migrate structured assemblers that only append unique entries.
3. Add row capacity reservation for `TreeLinearization`.
4. Remove avoidable temporary allocations in model coefficient evaluators.
5. Split static and dynamic coefficients so constant entries are not rewritten every Newton iteration.
6. Add a direct structured coefficient update path for serial `TreeNewtonLinearSolver`.
7. SIMD/batch dynamic model coefficient calculations if profiling still shows arithmetic cost.
8. Add validation tests comparing generic and direct structured assembly.

## Expected Outcome

The best near-term speedups should come from reducing generic assembly overhead:

- Less `set_value(...)` row searching.
- Fewer temporary vectors and heap allocations in model coefficient evaluators.
- Less repeated constant coefficient insertion.
- Less per-iteration coefficient resolve work.
- Better cache locality for coefficient values passed to the tree solver.

The SIMD/batch dense tree solve is already effective for the gen16 rigid case. Further speedup should focus on getting structured coefficients into the solver faster.
