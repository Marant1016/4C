# Last Code Changes

Date: 2026-08-17

## Goal

Support root-inlet flow boundary conditions in the optimized serial `NewtonTree` correction solve.

The tree solver previously found a boundary condition at the root inlet but always used that row to
compute the root inlet pressure correction directly. That worked for root pressure boundary
conditions, but a root flow boundary constrains the root inlet flow dof instead. The NOX and
`NewtonSparse` sparse workflows already assemble flow boundary conditions as ordinary sparse system
rows, so this change makes the tree root closure handle the same boundary type without adding any
sparse fallback to the production tree path.

## Code Updated

Updated:

```text
src/reduced_lung/src/solver/tree/4C_reduced_lung_tree_linear_solver.hpp
src/reduced_lung/src/solver/tree/4C_reduced_lung_tree_linear_solver.cpp
```

The tree solver now stores the actual root boundary type and constrained local dof from
`TreeBoundaryConditionMetadata`.

Root closure behavior is now:

```text
Pressure root BC:
  delta_p_root = rhs_root_boundary / coeff_root_boundary

Flow root BC:
  root bottom-up relation: delta_q_root = G_root * delta_p_root + h_root
  boundary equation:       coeff_root_boundary * delta_q_root = rhs_root_boundary
  delta_p_root = (rhs_root_boundary / coeff_root_boundary - h_root) / G_root
```

The flow case reuses the already computed bottom-up affine root relation, so the optimized
bottom-up/top-down algorithm and direct structured coefficient path remain unchanged. If the root
flow relation has a near-zero pressure slope, `NewtonTree` now fails fast with a specific error
instead of producing an undefined root pressure correction.

## Tests Updated

Updated:

```text
src/reduced_lung/tests/4C_reduced_lung_tree_linear_solver_test.cpp
```

Added a serial airway fixture with a root inlet flow boundary and outlet pressure boundary. New
coverage checks:

```text
TreeNewtonLinearSolver root-inlet-flow correction against SparseNewtonLinearSolver
direct structured coefficient assembly for root-inlet-flow trees against the generic structured path
full NewtonTree workflow with root inlet flow against NOX and NewtonSparse
```

## Validation Run

Passed:

```text
cmake --build build/debug --target reduced_lung_objs --parallel 4
cmake --build build/debug --target unittests_reduced_lung --parallel 4
ctest -R "^unittests_reduced_lung$" --output-on-failure
```

Date: 2026-08-09

## Goal

Make the reduced-lung `NewtonTree` workflow distributed-free and keep only the optimized serial tree solver path.

The distributed tree solver was an experimental MPI implementation that was not optimized to the same level as the serial tree solver. Keeping it in the final code would make the solver layer harder to read and could invite users to benchmark or rely on an unsupported path.

## Code Removed

Removed the distributed tree linear solver implementation:

```text
src/reduced_lung/src/solver/tree/4C_reduced_lung_distributed_tree_linear_solver.cpp
src/reduced_lung/src/solver/tree/4C_reduced_lung_distributed_tree_linear_solver.hpp
```

Removed the distributed tree linear solver NP2 test:

```text
src/reduced_lung/tests/solver/tree/4C_reduced_lung_distributed_tree_linear_solver_test.np2.cpp
```

Removed the dedicated distributed NewtonTree regression input:

```text
tests/input_files/reduced_lung_distributed_bifurcation_newton_tree.4C.yaml
```

Removed multi-rank NewtonTree regression registrations from:

```text
tests/list_of_tests.cmake
```

The serial `reduced_lung_lung_tree_gen10_newton_tree.4C.yaml` regression remains registered. Its removed `NP 2` registration is intentional because `NewtonTree` is no longer an MPI solver path.

The reduced-lung Google Benchmark folder had already been removed before this cleanup, so there was no active distributed benchmark source left to delete.

## Runtime Behavior

`NewtonTree` is now explicitly serial-only.

If a user runs a reduced-lung input with:

```yaml
nonlinear_solver: NewtonTree
```

using more than one MPI rank, setup fails immediately with a clear error message. The intended user action is to either:

```text
run NewtonTree with one MPI rank
```

or explicitly select a solver workflow that supports MPI runs:

```yaml
nonlinear_solver: NewtonSparse
```

or:

```yaml
nonlinear_solver: Nox
```

This fail-fast behavior is preferred over silently falling back to another solver because an automatic fallback would make benchmark results misleading.

## Code Updated

Updated:

```text
src/reduced_lung/src/4C_reduced_lung_main.cpp
```

The old runtime branch:

```text
comm_size == 1 -> TreeNewtonLinearSolver
comm_size > 1  -> DistributedTreeNewtonLinearSolver
```

was replaced by:

```text
comm_size != 1 -> throw a serial-only NewtonTree error
comm_size == 1 -> TreeNewtonLinearSolver
```

Updated:

```text
src/reduced_lung/src/solver/4C_reduced_lung_solver_profiles.hpp
```

Removed distributed-only tree profile counters:

```text
communication_time
communicated_coefficient_count
communicated_residual_count
boundary_relation_message_count
boundary_pressure_message_count
correction_scatter_message_count
communication_bytes
```

The remaining `TreeNewtonLinearSolverProfile` fields describe the serial tree solver and its optimized dense/SIMD/batch paths.

## Current Solver Policy

The final reduced-lung nonlinear solver options are still:

```text
Nox
NewtonSparse
NewtonTree
```

Their intended roles are:

```text
Nox          default/reference workflow, MPI-capable through existing infrastructure
NewtonSparse custom Newton workflow with sparse linear solver backend, MPI-capable
NewtonTree   optimized custom Newton workflow with serial tree linear solver backend only
```

`Nox` remains the safe default. `NewtonTree` remains an explicit opt-in path for serial performance runs.

## Relocated Test Inputs

The reduced-lung input files that were created for the custom Newton workflow but are not part of
the original upstream `tests/input_files/` set have been moved out of the 4C test-input directory.

They are now kept as a clean, self-contained archive in:

```text
/scratch/Rodriguez/workspace/4C/files/unit_tests/
```

This folder includes the moved `NewtonSparse`, `NewtonTree`, and generated gen10 reduced-lung YAML
inputs, plus the supporting `*_fields.json` files needed by their relative `from_file` references.

The corresponding `four_c_test(...)` registrations were removed from:

```text
tests/list_of_tests.cmake
```

so CTest no longer expects those files under `tests/input_files/`.

## Relocated GoogleTests and Benchmarks

The reduced-lung GoogleTest source files created during the custom Newton/NewtonTree work were moved
out of the auto-discovered 4C test tree and archived in:

```text
/scratch/Rodriguez/workspace/4C/files/google_tests/reduced_lung/
```

This includes the Newton solver, Newton-vs-NOX, tree metadata, tree linear solver, and input-pipeline
test sources. The distributed `.np2` Newton test sources were archived there as well because
`NewtonTree` is now serial-only and distributed NewtonTree tests would be misleading.

## Restored Solver GoogleTests

The main custom Newton and serial tree-solver GoogleTests have been restored into the active
auto-discovered reduced-lung test folder:

```text
src/reduced_lung/tests/4C_reduced_lung_newton_solver_test.np2.cpp
src/reduced_lung/tests/4C_reduced_lung_tree_metadata_test.cpp
src/reduced_lung/tests/4C_reduced_lung_tree_linear_solver_test.cpp
```

The restored coverage verifies:

```text
NewtonSparse custom Newton workflow on an analytical terminal-unit problem
tree metadata construction and validation invariants
TreeNewtonLinearSolver corrections against the sparse Newton linear solver
structured tree coefficient assembly against the generic structured path
small NewtonTree workflows against NOX and NewtonSparse reference workflows
```

These tests intentionally live in the flat `src/reduced_lung/tests/` directory to match the current
`four_c_auto_define_tests()` layout used by the existing reduced-lung unit tests.

## Restored Runtime Input Smoke Tests

Four small reduced-lung runtime input tests have been restored into the active input-file test set:

```text
tests/input_files/reduced_lung_terminal_unit_newton_sparse.4C.yaml
tests/input_files/reduced_lung_terminal_unit_newton_tree.4C.yaml
tests/input_files/reduced_lung_aw_bifurcation_flow_newton_sparse.4C.yaml
tests/input_files/reduced_lung_aw_bifurcation_flow_newton_tree.4C.yaml
```

They reuse the existing supporting field files:

```text
tests/input_files/reduced_lung_terminal_unit_fields.json
tests/input_files/reduced_lung_aw_bifurcation_fields.json
```

The tests are registered in:

```text
tests/list_of_tests.cmake
```

These are intentionally lightweight serial smoke tests. Their role is to verify that real `.4C.yaml`
inputs select and execute `nonlinear_solver: NewtonSparse` and `nonlinear_solver: NewtonTree` through
the normal 4C runtime input pipeline. Deeper solver correctness remains covered by the restored
GoogleTests in `src/reduced_lung/tests/`.

The generated gen16 benchmark inputs are kept in:

```text
/scratch/Rodriguez/workspace/4C/files/gen16_inputs/
```

The benchmark scripts and historical benchmark outputs are kept in:

```text
/scratch/Rodriguez/workspace/4C/files/benchmarks/scripts/
/scratch/Rodriguez/workspace/4C/files/benchmarks/results/
```

## Validation To Run

After this cleanup, run:

```text
cmake --build build/debug --target reduced_lung_objs --parallel 4
cmake --build build/debug --target unittests_reduced_lung unittests_reduced_lung.np2 --parallel 4
ctest -R "^unittests_reduced_lung(\.np2)?$" --output-on-failure
ctest -R "^reduced_lung_.*newton_.*\.4C\.yaml-p1$" --output-on-failure
```

For large end-to-end benchmarking, the currently used generated-tree scripts remain independent of the removed Google Benchmark sources:

```text
RUNS=3 /scratch/Rodriguez/workspace/4C/files/benchmarks/scripts/benchmark_reduced_lung_gen16_linear_newton_tree.sh
RUNS=3 /scratch/Rodriguez/workspace/4C/files/benchmarks/scripts/benchmark_reduced_lung_gen16_nonlinear_newton_tree.sh
RUNS=3 /scratch/Rodriguez/workspace/4C/files/benchmarks/scripts/benchmark_reduced_lung_gen16_solvers.sh
```
