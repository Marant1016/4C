# Assemply Step 1

## Goal

Add fine-grained timing for structured tree linearization assembly. This step does not optimize assembly behavior yet. It only splits the existing `tree_assembly_s` profile time into actionable sub-counters so the next optimization can target the real bottleneck.

## Motivation

The gen16 YAML profile showed:

```text
tree_assembly_s: 15.9798
tree_solve_s: 4.17204
tree_dense_s: 0.687825
tree_simd_lanes: 92459000
tree_dense_fallbacks: 0
tree_unsupported_fallbacks: 0
```

This proves the SIMD/batch tree solve is already active and fallback-free. The largest measured cost is now structured tree linearization assembly, so this step instruments that phase in more detail.

## Code Changes

### Named Tree Assembly Phases

File:

```text
src/reduced_lung/src/4C_reduced_lung_helpers.hpp
```

Added `ReducedLungAssemblyPipeline::TreeLinearizationAssemblyPhase`:

```cpp
enum class TreeLinearizationAssemblyPhase
{
  Airways,
  TerminalUnits,
  Junctions,
  BoundaryConditions,
  Other,
};
```

Added `NamedTreeLinearizationAssembler` so each tree-linearization callback carries a phase label:

```cpp
struct NamedTreeLinearizationAssembler
{
  TreeLinearizationAssemblyPhase phase = TreeLinearizationAssemblyPhase::Other;
  TreeLinearizationAssembler callback;
};
```

The pipeline now stores:

```cpp
std::vector<NamedTreeLinearizationAssembler> tree_linearization_assemblers;
```

instead of a vector of anonymous callbacks.

### Pipeline Registration

File:

```text
src/reduced_lung/src/4C_reduced_lung_helpers.cpp
```

The four structured tree linearization callbacks are now registered with explicit phase labels:

```text
Airways
TerminalUnits
Junctions
BoundaryConditions
```

The callback behavior is unchanged. Only the metadata around each callback changed.

### New Profile Counters

File:

```text
src/reduced_lung/src/4C_reduced_lung_solver_profile.hpp
```

Added these fields to `NewtonSolverProfile`:

```cpp
double tree_linearization_clear_time = 0.0;
double tree_linearization_airway_time = 0.0;
double tree_linearization_terminal_unit_time = 0.0;
double tree_linearization_junction_time = 0.0;
double tree_linearization_boundary_condition_time = 0.0;
double tree_linearization_other_time = 0.0;
double tree_linearization_solver_update_time = 0.0;
```

The existing total remains unchanged:

```cpp
double structured_tree_linearization_assembly_time = 0.0;
```

### Timing Implementation

File:

```text
src/reduced_lung/src/4C_reduced_lung_newton_solver.cpp
```

Updated `NewtonSolver::assemble_tree_linearization_for_current_state()` to time:

```text
TreeLinearization clear/reset
Airway callback
Terminal-unit callback
Junction callback
Boundary-condition callback
Other callback, if any
linear_solver_->set_tree_linearization(...)
```

The total `structured_tree_linearization_assembly_time` still wraps the whole function, so previous reporting remains comparable.

### Real YAML Profile Output

File:

```text
src/reduced_lung/src/4C_reduced_lung_main.cpp
```

The `FOUR_C_REDUCED_LUNG_TREE_PROFILE=1` summary now prints:

```text
tree_assembly_clear_s
tree_assembly_airways_s
tree_assembly_terminal_units_s
tree_assembly_junctions_s
tree_assembly_boundary_conditions_s
tree_assembly_other_s
tree_assembly_solver_update_s
```

alongside the existing total:

```text
tree_assembly_s
```

### Benchmark Counter Output

Files:

```text
src/reduced_lung/benchmark_tests/4C_reduced_lung_solver_benchmark.cpp
src/reduced_lung/benchmark_tests/4C_reduced_lung_distributed_tree_solver_benchmark.cpp
```

The same fine-grained assembly counters are exposed in benchmark output.

### Call Site Updates

Updated all existing direct loops over `tree_linearization_assemblers` to call:

```cpp
tree_linearization_assembler.callback(...)
```

Affected files:

```text
src/reduced_lung/benchmark_tests/4C_reduced_lung_solver_benchmark.cpp
src/reduced_lung/benchmark_tests/4C_reduced_lung_distributed_tree_solver_benchmark.cpp
src/reduced_lung/tests/4C_reduced_lung_tree_linear_solver_test.cpp
src/reduced_lung/tests/4C_reduced_lung_distributed_tree_linear_solver_test.np2.cpp
```

## How To Use The New Counters

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

Look at the final profile block. The new counters identify where `tree_assembly_s` is spent.

Interpretation guide:

- High `tree_assembly_clear_s`: optimize `TreeLinearization::clear_values()` and row storage reuse.
- High `tree_assembly_airways_s`: optimize airway derivative computation, temporary allocations, and coefficient writes.
- High `tree_assembly_terminal_units_s`: optimize terminal-unit gradient/rheology coefficient loops.
- High `tree_assembly_junctions_s`: optimize junction coefficient insertion and static coefficient reuse.
- High `tree_assembly_boundary_conditions_s`: optimize boundary coefficient insertion/static reuse.
- High `tree_assembly_solver_update_s`: optimize `TreeNewtonLinearSolver::set_tree_linearization(...)` and `resolve_structured_coefficient_locations(...)`.

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

Run the gen16 YAML profile again and inspect the new breakdown. The next optimization should be selected from the largest sub-counter, not from the old aggregate `tree_assembly_s` alone.
