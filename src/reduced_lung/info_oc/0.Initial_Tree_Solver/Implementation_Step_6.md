# Implementation Step 6

Date: 2026-07-18

## Scope

Implemented a standalone reduced-lung tree metadata builder for future tree-based Newton linear solvers.

This step only builds and validates metadata. It does not implement the tree linear solver, does not compute bottom-up/top-down Schur-complement relations, and is not wired into the production NOX runtime path or the custom Newton path yet.

The implementation uses the Jupyter notebook `05_airway_tree_solver_kv_soa_numpy_nonlinear.ipynb` as conceptual guidance for preprocessing the tree once into parent/child relations and traversal layers. The C++ metadata is adapted to the current 4C reduced-lung formulation, which uses duplicated endpoint dofs, explicit junction equations, explicit boundary-condition equations, and existing global/local row and dof maps.

The `src/reduced_lung/src/1d_pipe_flow/` implementation was intentionally ignored.

## Files Added

### `src/reduced_lung/src/4C_reduced_lung_tree_metadata.hpp`

Added public metadata types for the future tree solver:

- `TreeElementKind`
- `TreeJunctionKind`
- `TreeBoundarySide`
- `TreeElementMetadata`
- `TreeJunctionMetadata`
- `TreeBoundaryConditionMetadata`
- `ReducedLungTreeMetadata`
- `ReducedLungTreeMetadataContext`

Added the builder entry point:

```text
build_reduced_lung_tree_metadata(context)
```

`ReducedLungTreeMetadataContext` collects metadata inputs from the existing setup:

- `ReducedLungParameters`, including input topology,
- `first_global_dof_of_ele`,
- `global_dof_per_ele`,
- airway model containers,
- terminal-unit model containers,
- connection and bifurcation containers,
- boundary-condition containers,
- row map,
- locally relevant dof map.

### `src/reduced_lung/src/4C_reduced_lung_tree_metadata.cpp`

Implemented the metadata builder.

For each element, the builder records:

- global element id,
- element kind, airway or terminal unit,
- inlet and outlet node ids,
- parent element index,
- up to two child element indices,
- first global dof id,
- number of dofs,
- global dof ids,
- local dof ids on the locally relevant dof map,
- first local/global state-equation ids,
- number of state equations.

For each junction, the builder records:

- connection or bifurcation kind,
- parent element index,
- child element indices,
- first local/global junction-equation ids,
- number of junction equations,
- global and local dof ids used by the junction equations.

For each boundary condition, the builder records:

- boundary type,
- inlet or outlet side,
- attached node id,
- attached element index,
- local/global equation id,
- global/local constrained dof id.

The builder also records:

- `root_element_index`,
- `root_node_id`,
- airway element indices,
- terminal-unit element indices,
- bottom-up traversal layers,
- top-down traversal layers,
- global dof count,
- global equation count,
- locally relevant dof count.

## Validation Added

The builder validates the topology and metadata before returning.

Current checks:

- topology has at least one element,
- each topology element has exactly two nodes,
- topology node ids are valid 1-based input ids,
- no self-connected element exists,
- each element has a known dof offset and dof count,
- each element has model equation metadata,
- each element has at least three dofs,
- global equation count equals global dof count,
- each directed child has at most one parent,
- branch degree is at most two children,
- directed topology is acyclic,
- exactly one root element exists,
- all elements are reachable from the root,
- terminal units are leaves,
- connection metadata matches the directed topology,
- bifurcation metadata matches the directed topology,
- every parent with one child has a connection,
- every parent with two children has a bifurcation,
- every boundary condition is attached to the inlet or outlet side of its element,
- every boundary condition constrains the expected dof for its side and type,
- the root inlet side has a boundary condition,
- every leaf outlet side has a boundary condition.

## Notebook Adaptation

The notebook prototype stores solver work arrays such as:

- `G`,
- `h`,
- pressure corrections,
- flow corrections,
- residuals,
- tangents.

Those arrays were intentionally not added in Phase 6. They belong to Phase 7, when the actual tree linear solver is implemented.

Phase 6 only stores stable structural metadata that can be reused by a future solver implementation.

The notebook builds layers partly from generation data. The C++ implementation instead builds traversal layers from directed connectivity:

```text
element_nodes[0] = inlet / parent side
element_nodes[1] = outlet / child side
```

This avoids trusting generation metadata for solver traversal.

## Tests Added

Added:

```text
src/reduced_lung/tests/4C_reduced_lung_tree_metadata_test.cpp
```

The new tests cover:

- valid connection metadata,
- valid bifurcation metadata with terminal-unit leaves,
- bottom-up and top-down layer construction,
- airway and terminal-unit element classification,
- junction dof metadata,
- boundary-side classification,
- directed cycle validation,
- unsupported branch degree validation,
- missing root boundary validation.

The tests construct small metadata fixtures directly from reduced-lung containers and maps. They do not run NOX, custom Newton, residual assembly, or Jacobian assembly.

## Preserved Behavior

- The production runtime path still uses `NoxSolver`.
- The custom Newton path is unchanged.
- No YAML input behavior was changed.
- No residual equation was changed.
- No Jacobian equation was changed.
- No row ordering was changed.
- No dof ordering was changed.
- No tree linear solver was added yet.

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

- Metadata is not yet consumed by a solver.
- The tree-based linear solver is still not implemented.
- No tree-solver work arrays are allocated yet.
- No solver selection option was added to YAML input.
- The first metadata implementation is designed for the current duplicated endpoint dof formulation.
