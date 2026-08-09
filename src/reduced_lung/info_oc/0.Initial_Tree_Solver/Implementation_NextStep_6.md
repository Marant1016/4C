# Implementation NextStep 6 - Subtree-Partition Communication

## Scope

This step replaces the Step 5 replicated distributed tree solve with rank-local tree work and small
partition-boundary messages.

The serial path is unchanged:

- `comm_size == 1` still uses `TreeNewtonLinearSolver`.
- `comm_size > 1` still uses `DistributedTreeNewtonLinearSolver`.
- `NewtonSolver`, `NewtonSparse`, and `Nox` remain unchanged.

## Files Modified

- `src/reduced_lung/src/4C_reduced_lung_distributed_tree_linear_solver.cpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_solver_profile.hpp`
- `src/reduced_lung/tests/4C_reduced_lung_distributed_tree_linear_solver_test.np2.cpp`

## Distributed Algorithm

The distributed tree solver no longer gathers all residual entries or all structured coefficients.
Instead, each rank solves only the elements it owns according to `TreeElementMetadata::owner_rank`.

The solver still uses the same Newton correction convention:

```text
J * delta = -F
x_new = x_old + delta
```

### Local Coefficient Access

The solver now uses a rank-local coefficient provider:

```text
global row id -> local residual row id
global dof id -> local structured-tree column id
```

Owned element blocks are assembled only from locally available rows and locally relevant dofs. If a
layout does not provide a row/dof required by an owned element block, the solver fails fast with a
clear error.

### Bottom-Up Phase

For each bottom-up tree layer:

1. A rank processes only elements it owns.
2. Same-rank child subtree relations are read directly from local storage.
3. Cross-rank child subtree relations are received as boundary messages.
4. The rank condenses each owned element to:

```text
q_in = G * p_in + h
```

5. If the parent element is owned by another rank, the child owner sends:

```text
element_index, G, h
```

These messages are exchanged with `MPI_Alltoallv` per tree layer. Only partition-boundary relations
are sent.

### Top-Down Phase

The root owner computes the root inlet-pressure correction from the root inlet boundary row.

For each top-down tree layer:

1. A rank processes only owned elements whose inlet pressure is known.
2. The rank recovers local element correction values.
3. Same-rank child inlet pressures are stored locally.
4. Cross-rank child inlet pressures are sent as:

```text
child_element_index, child_inlet_pressure
```

These messages are also exchanged with `MPI_Alltoallv` per tree layer.

### Correction Scatter

Element owners compute correction values, but the Newton correction vector is distributed by the
correction map. Therefore, the solver now performs a final correction scatter:

```text
global_dof_id, delta_value
```

If the computed dof belongs to the local correction map, it is inserted directly. Otherwise it is sent
to the correction-vector owner and inserted there.

## Removed Step 5 Bottleneck

The distributed tree solver no longer performs:

- global residual gather
- global structured-coefficient gather
- replicated full-tree solve on every rank

Communication is now limited to:

- boundary subtree relations `(G, h)`
- boundary inlet-pressure corrections
- final correction-vector scatter messages

## Profiling Additions

`TreeNewtonLinearSolverProfile` now tracks the new message categories:

- `boundary_relation_message_count`
- `boundary_pressure_message_count`
- `correction_scatter_message_count`

The existing `communication_time` and `communication_bytes` counters are updated from the new
boundary exchanges.

## Validation Added

The distributed tree np2 test was strengthened:

- Serial airway chain correction compared against `SparseNewtonLinearSolver`.
- Bifurcation airway correction compared against `SparseNewtonLinearSolver`.
- Both cases assert that the generated tree metadata actually contains at least one cross-rank edge.

Relevant tests:

```text
ReducedLungDistributedTreeLinearSolverTests.SerialAirwaysMatchSparseSolverOnTwoRanks
ReducedLungDistributedTreeLinearSolverTests.BifurcationAirwaysMatchSparseSolverOnTwoRanks
```

## Verification Run

The following commands were run successfully during this step:

```bash
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake \
  --build build/debug \
  --target unittests_reduced_lung \
  --parallel 4
```

```bash
/scratch/Rodriguez/workspace/CLion-2026.1.3/clion-2026.1.3/bin/cmake/linux/x64/bin/cmake \
  --build build/debug \
  --target unittests_reduced_lung.np2 \
  --parallel 4
```

```bash
ctest -R "^unittests_reduced_lung$" --output-on-failure
```

```bash
ctest -R "^unittests_reduced_lung\.np2$" --output-on-failure
```

## Current Limitations

- The solver still uses collective `MPI_Alltoallv` exchanges per tree layer rather than fully
  asynchronous point-to-point sends and receives.
- The correction owner table is queried during each solve from the correction map.
- The implementation still assumes that each owned element's rows and locally relevant dofs are
  available on the element owner rank.
- More complex distributed mixed airway/terminal-unit layouts should be added to np2 validation next.

## Natural Next Step

Profile the distributed tree path on larger np2/np4 cases and reduce synchronization overhead:

- cache correction owner lookup in the distributed solver
- precompute per-layer send/receive plans
- replace per-layer all-to-all exchanges with point-to-point communication between actual neighboring
  tree partitions
- add mixed airway/terminal-unit distributed correction and full-workflow tests
