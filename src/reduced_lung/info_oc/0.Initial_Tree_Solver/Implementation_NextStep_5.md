# Implementation NextStep 5 - Parallelization Strategy

## Scope

This step adds the first MPI-capable `NewtonTree` path without changing the existing nonlinear
Newton workflow or the validated serial structured-block tree solver.

The implementation is intentionally conservative:

- `comm_size == 1` keeps using the existing `TreeNewtonLinearSolver`.
- `comm_size > 1` uses a new `DistributedTreeNewtonLinearSolver`.
- The distributed solver gathers locally assembled structured tree coefficients and residual values,
  performs the tree solve in global IDs, and writes back only locally owned correction entries.
- `NewtonSparse` and `Nox` remain unchanged reference/fallback workflows.

This is not yet subtree-partitioned bottom-up/top-down communication. It is the foundation that makes
the structured tree path work with distributed reduced-lung maps while preserving the serial
algorithm's mathematics.

## Files Added

- `src/reduced_lung/src/4C_reduced_lung_distributed_tree_linear_solver.cpp`
- `src/reduced_lung/tests/4C_reduced_lung_distributed_tree_linear_solver_test.np2.cpp`
- `src/reduced_lung/info_oc/Implementation_NextStep_5.md`

## Files Modified

- `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_metadata.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_metadata.cpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linearization.hpp`
- `src/reduced_lung/src/4C_reduced_lung_tree_linearization.cpp`
- `src/reduced_lung/src/4C_reduced_lung_solver_profile.hpp`
- `src/reduced_lung/src/4C_reduced_lung_main.cpp`

## Main Changes

### Distributed Tree Metadata

`ReducedLungTreeMetadata` now stores owner ranks for:

- tree elements
- junction equations
- boundary-condition equations

Metadata construction now all-reduces local element equation metadata, junction metadata, and boundary
metadata so every rank can build the same global tree topology. Local row/dof IDs remain valid on the
owning/relevant rank and may be `-1` elsewhere.

### Distributed Tree Linear Solver

Added `DistributedTreeNewtonLinearSolver` and `DistributedTreeNewtonLinearSolverContext`.

The solver uses the same correction convention as the serial path:

```text
J * delta = -F
x_new = x_old + delta
```

The current MPI implementation does this per Newton correction:

1. Gather distributed residual entries into a global residual vector.
2. Gather local structured tree-linearization coefficients as `(global_row, global_dof, value)`.
3. Build a global coefficient provider.
4. Run the existing tree condensation/recovery algorithm in global IDs.
5. Insert only correction entries present in the local `delta` map.

This keeps assembly ownership unchanged and avoids sparse Jacobian assembly for distributed
`NewtonTree` runs.

### Runtime Selection

`ReducedLungSimulation::build_newton_solver_with_tree_linear_solver()` now selects:

- `TreeNewtonLinearSolver` for serial runs.
- `DistributedTreeNewtonLinearSolver` for MPI runs.

The previous hard error for `NewtonTree` with `comm_size > 1` was removed.

### Tree Linearization Access

`TreeLinearization` now exposes row entries through:

```cpp
const std::vector<std::pair<int, double>>& entries(int local_row_id) const;
```

This is used by the distributed solver to gather structured coefficients without converting through a
sparse matrix.

### Profiling

`TreeNewtonLinearSolverProfile` now includes distributed communication counters:

- `communication_time`
- `communicated_coefficient_count`
- `communicated_residual_count`
- `communication_bytes`

The new distributed solver records these counters when a profile is attached.

## Validation Added

Added a two-rank unit test:

```text
ReducedLungDistributedTreeLinearSolverTests.SerialAirwaysMatchSparseSolverOnTwoRanks
```

The test assembles a distributed serial-airway reduced-lung system, solves one Newton correction with
`SparseNewtonLinearSolver`, solves the same correction with `DistributedTreeNewtonLinearSolver`, and
compares the distributed correction vectors.

## Verification Run

The following commands were run successfully:

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

Both unit-test targets passed.

## Current Limitations

- The distributed solver currently gathers and replicates the complete structured correction solve on
  all ranks.
- It does not yet partition subtrees and exchange only `(G, h)` and inlet-pressure messages across
  cut edges.
- This first MPI path prioritizes correctness and compatibility with existing assembly maps over
  distributed scalability.

## Natural Next Step

Replace the replicated global solve with true partition-boundary communication:

- local bottom-up condensation per rank
- send `(G, h)` from child partition to parent partition
- solve the root relation
- send inlet-pressure corrections from parent partition to child partition
- recover local corrections without gathering all coefficients globally
