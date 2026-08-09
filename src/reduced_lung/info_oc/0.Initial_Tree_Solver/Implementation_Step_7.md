# Implementation Step 7

Date: 2026-07-18

## Scope

Implemented the first serial tree-based Newton linear solver for the reduced-lung model.

This solver is exposed through the Phase 5 `NewtonLinearSolver` interface and consumes the Phase 6 `ReducedLungTreeMetadata`. It uses the already assembled sparse Jacobian and residual, so existing residual and Jacobian physics remain unchanged.

The implementation performs bottom-up subtree condensation and top-down correction recovery on the current duplicated endpoint dof formulation.

The production runtime path still uses NOX. The custom Newton path is not selectable from YAML yet. The new tree solver is tested directly and is not wired into `ReducedLungSimulation`.

The `src/reduced_lung/src/1d_pipe_flow/` implementation was intentionally ignored.

## Files Added

### `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp`

Added:

- `TreeNewtonLinearSolverContext`
- `TreeNewtonLinearSolver`

`TreeNewtonLinearSolver` implements:

```text
NewtonLinearSolver::solve(jacobian, residual, x, metadata, delta)
```

The sign convention remains the Phase 5 convention:

```text
J * delta = -F
x_new = x_old + delta
```

The context currently stores:

- a reference to `ReducedLungTreeMetadata`,
- a dense local pivot tolerance.

### `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp`

Implemented the serial tree linear solver.

Current restrictions:

- serial only,
- requires a completed sparse Jacobian,
- requires all residual rows and correction dofs to be locally available,
- requires Phase 6 tree metadata,
- requires one outlet boundary for each leaf element,
- requires one root inlet boundary,
- supports up to two children per parent, matching connection and bifurcation metadata.

## Algorithm

The implementation follows the same conceptual structure as the Jupyter notebook, adapted to the 4C formulation.

Notebook concept:

```text
bottom-up: compute dQ = G dPin + h
top-down: recover dPin, dPout, dQ
```

C++ implementation:

```text
bottom-up: compute one subtree relation per element
top-down: recover all global dof corrections into delta
```

For every element, the solver computes a scalar inlet relation:

```text
delta_q_in = G * delta_p_in + h
```

The relation is assembled from existing sparse Jacobian coefficients and the RHS `-residual`.

### Bottom-Up Condensation

Elements are processed using `ReducedLungTreeMetadata::bottom_up_layers`.

For each element, the solver builds a small dense local system with all element dofs except inlet pressure as local unknowns. The inlet pressure correction is treated as the remaining subtree parameter.

Leaf element local rows:

- element state-equation rows,
- the leaf outlet boundary-condition row.

Internal element local rows:

- element state-equation rows,
- the junction flow-conservation row.

Pressure-continuity junction rows are used to substitute each child inlet pressure in terms of the parent outlet pressure. Each child inlet flow is then substituted using the child subtree relation.

The small dense system is solved twice:

- once for the constant/intercept part,
- once for the slope with respect to inlet pressure.

The inlet-flow entry gives the element subtree relation `(G, h)`.

### Top-Down Recovery

The root inlet boundary condition gives the root inlet pressure correction.

Elements are then processed using `ReducedLungTreeMetadata::top_down_layers`.

For each element, the stored local affine recovery data reconstructs all element dof corrections, which are written into the global `delta` vector using the existing global dof ids.

For each child, the pressure-continuity row gives the child inlet pressure correction from the recovered parent outlet pressure correction.

### Dense Local Solver

Added a small local Gaussian-elimination helper with partial pivoting inside the tree solver implementation file.

It is intentionally local to this first tree solver implementation because all condensed blocks are small.

The solver throws if a local pivot is below `pivot_tolerance`, which catches singular or underconstrained local subtree blocks.

## Notebook Adaptation

The notebook stores model-specific arrays such as:

- `aw_G`,
- `aw_h`,
- `tu_G`,
- `tu_h`,
- pressure corrections,
- flow corrections.

The C++ implementation keeps only generic solver work arrays:

- one `SubtreeRelation` per element,
- one `ElementRecoveryData` block per element,
- one inlet-pressure correction per element during top-down recovery.

This avoids hard-coding the notebook's simplified Kelvin-Voigt-only data layout. The first C++ solver instead reads coefficients from the assembled Jacobian, so it can work with the currently assembled reduced-lung element, junction, and boundary equations as long as they fit the validated tree structure.

## Tests Added

Added:

```text
src/reduced_lung/tests/4C_reduced_lung_tree_linear_solver_test.cpp
```

The new tests assemble reduced-lung residuals and Jacobians through the existing assembly pipeline, then compare one Newton correction solve between:

- `SparseNewtonLinearSolver`, backed by `Core::LinAlg::Solver` / UMFPACK,
- `TreeNewtonLinearSolver`, backed by the new serial bottom-up/top-down tree algorithm.

Test cases added:

- single terminal unit,
- serial rigid airways,
- rigid-airway bifurcation.

The tests compare the full correction vector `delta` entry-by-entry.

## Preserved Behavior

- The existing NOX runtime path was not changed.
- The custom Newton solver interface was not changed.
- No YAML input behavior was changed.
- No residual equation was changed.
- No Jacobian equation was changed.
- No dof ordering was changed.
- No row ordering was changed.
- No Phase 6 tree metadata validation was relaxed.

## Verification Performed

Built the serial reduced-lung unit-test executable:

```text
cmake --build build/debug --target unittests_reduced_lung --parallel 4
```

Result: passed.

Ran the serial reduced-lung unit-test target:

```text
ctest -R "^unittests_reduced_lung$" --output-on-failure
```

Result: passed, 1 of 1 selected test target.

Ran the serial and `.np2` reduced-lung unit-test targets:

```text
ctest -R "^unittests_reduced_lung(\.np2)?$" --output-on-failure
```

Result: passed, 2 of 2 selected test targets.

## Current Limitations

- The tree solver is serial-only.
- The tree solver is not selectable from YAML input.
- The production runtime still uses NOX.
- The tree solver currently consumes the assembled sparse Jacobian instead of structured model-derivative blocks.
- Solver diagnostics are limited to clear exceptions on unsupported parallel execution, missing closure, missing/near-zero required coefficients, and singular local dense pivots.
- Parallel tree elimination is not implemented.
