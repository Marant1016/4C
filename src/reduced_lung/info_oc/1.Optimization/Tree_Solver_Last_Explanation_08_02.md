# Tree Solver Current Explanation - 2026-08-02

This document explains the current tree-based Newton linear solver in `src/reduced_lung/` as of 2026-08-02. It focuses on the tree solver and its integration into the reduced-lung custom Newton workflow, not on the full 4C framework.

The most important point is:

```text
NewtonTree = custom full-step Newton loop + tree-structured Newton correction solve
```

The tree solver does not change the reduced-lung equations. It solves the same linearized Newton correction system as the sparse solver,

```text
J * delta = -F
x_new = x_old + delta
```

but exploits the directed airway-tree topology instead of factoring the full sparse matrix.

## Overview

The reduced-lung model uses a duplicated-endpoint formulation. Each element owns its own pressure and flow degrees of freedom. Neighboring elements are not coupled by shared nodal unknowns. Instead, explicit junction rows enforce pressure continuity and flow conservation.

The tree solver uses this structure directly:

```text
leaf elements
  -> condensed into affine inlet relations
  -> parent elements absorb child relations
  -> repeat bottom-up until the root
  -> root boundary gives root inlet pressure correction
  -> corrections are recovered top-down
  -> global Newton correction vector delta is filled
```

For each subtree, the bottom-up pass computes the affine relation:

```text
delta_q_in = G * delta_p_in + h
```

where:

- `delta_p_in` is the inlet pressure correction of the subtree root element.
- `delta_q_in` is the inlet flow correction of the same element.
- `G` is the condensed slope of the subtree.
- `h` is the residual-dependent intercept.

This is a Schur-complement elimination over the directed tree. Children are eliminated before their parents, so each child subtree enters the parent solve only through its scalar `(G, h)` relation.

## Main Files and Classes

### Runtime Solver Selection

| File | Main role |
| --- | --- |
| `src/reduced_lung/src/4C_reduced_lung_main.cpp` | Builds reduced-lung simulation objects, selects `Nox`, `NewtonSparse`, or `NewtonTree`, and calls the selected solver each timestep. |
| `src/reduced_lung/src/4C_reduced_lung_input.hpp` | Defines `ReducedLungParameters::NonlinearSolverType` with values `Nox`, `NewtonSparse`, and `NewtonTree`. |
| `src/reduced_lung/src/4C_reduced_lung_input.cpp` | Registers the `nonlinear_solver` input parameter and keeps `Nox` as the default. |

### Shared Assembly and Newton Infrastructure

| File | Main role |
| --- | --- |
| `src/reduced_lung/src/4C_reduced_lung_helpers.hpp/cpp` | Defines `ReducedLungAssemblyPipeline` and registers residual, sparse Jacobian, structured tree-linearization, and state-update callbacks. |
| `src/reduced_lung/src/4C_reduced_lung_newton_solver.hpp/cpp` | Implements `NewtonSolver`, the custom full-step Newton loop used by `NewtonSparse` and `NewtonTree`. |
| `src/reduced_lung/src/4C_reduced_lung_linear_solver.hpp/cpp` | Defines `NewtonLinearSolver`, `NewtonLinearizationType`, `NewtonLinearSystemMetadata`, and `SparseNewtonLinearSolver`. |

### Tree Solver Infrastructure

| File | Main role |
| --- | --- |
| `src/reduced_lung/src/4C_reduced_lung_tree_metadata.hpp/cpp` | Builds and validates directed tree metadata from existing reduced-lung setup data. |
| `src/reduced_lung/src/4C_reduced_lung_tree_linearization.hpp/cpp` | Stores structured local derivative coefficients for `NewtonTree`. |
| `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp` | Declares the serial and distributed tree linear solvers. |
| `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp` | Implements the serial bottom-up/top-down tree solver. |
| `src/reduced_lung/src/4C_reduced_lung_distributed_tree_linear_solver.cpp` | Implements the MPI-capable tree solver using partition-boundary messages. |

### Physics Modules Supplying Residuals and Derivatives

| Component | Files |
| --- | --- |
| Airways | `src/reduced_lung/src/airways/4C_reduced_lung_airways.cpp`, `src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.cpp`, and related airway files. |
| Terminal units | `src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit.cpp`, `src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp`, and related terminal-unit files. |
| Junctions | `src/reduced_lung/src/4C_reduced_lung_junctions.hpp/cpp`. |
| Boundary conditions | `src/reduced_lung/src/4C_reduced_lung_boundary_conditions.hpp/cpp`. |

### Important Classes and Structs

| Type | Meaning |
| --- | --- |
| `ReducedLungAssemblyPipeline` | Ordered callback registry for residual assembly, sparse Jacobian assembly, structured tree-linearization assembly, and state updates. |
| `NewtonSolver` | Custom nonlinear Newton driver. |
| `NewtonLinearSolver` | Interface for one Newton correction solve. |
| `SparseNewtonLinearSolver` | Custom Newton sparse backend using `Core::LinAlg::Solver`. |
| `TreeNewtonLinearSolver` | Serial tree Newton correction solver. |
| `DistributedTreeNewtonLinearSolver` | MPI tree Newton correction solver. |
| `TreeLinearization` | Row-oriented local derivative storage for structured tree solves. |
| `TreeElementMetadata` | Per-element topology, dof, row, kind, and ownership metadata. |
| `TreeJunctionMetadata` | Per-junction connection or bifurcation metadata. |
| `TreeBoundaryConditionMetadata` | Per-boundary-condition row, side, dof, type, and ownership metadata. |
| `ReducedLungTreeMetadata` | Full immutable tree metadata consumed by tree solvers. |

## Solver Selection Workflow

Runtime selection is controlled by:

```yaml
REDUCED DIMENSIONAL LUNG DYNAMIC:
  nonlinear_solver: Nox        # default if omitted
```

The supported values are:

| Option | Nonlinear driver | Linear correction solver | Current purpose |
| --- | --- | --- | --- |
| `Nox` | Existing NOX wrapper | `Core::LinAlg::Solver` through NOX | Default/reference path. |
| `NewtonSparse` | `NewtonSolver` | `SparseNewtonLinearSolver` | Custom Newton sparse reference/debug path. |
| `NewtonTree` | `NewtonSolver` | `TreeNewtonLinearSolver` or `DistributedTreeNewtonLinearSolver` | Structured tree solver path. |

The switch happens in `ReducedLungSimulation::build_linear_system_and_solver()` in `4C_reduced_lung_main.cpp`.

The construction logic is:

```text
Nox
  -> build_nox_solver()

NewtonSparse
  -> build_newton_solver_with_sparse_linear_solver()
  -> SparseNewtonLinearSolver
  -> NewtonSolver

NewtonTree
  -> build_newton_solver_with_tree_linear_solver()
  -> build_reduced_lung_tree_metadata(...)
  -> if comm_size == 1:
       TreeNewtonLinearSolver with StructuredTreeBlocks
     else:
       DistributedTreeNewtonLinearSolver
  -> NewtonSolver
```

`Nox` remains the default. `NewtonSparse` and `NewtonTree` are explicit opt-in workflows.

Even in `NewtonTree`, the common sparse matrix object `sysmat_` is still allocated because it is part of the shared linear-solver interface. However, runtime `NewtonTree` does not assemble or complete the sparse Jacobian because the tree linear solver reports that it needs `StructuredTreeBlocks`.

## Newton Loop Integration

The custom Newton loop is implemented by `NewtonSolver::solve(double time)`.

For each physical timestep, it performs:

```text
current_time = time
for iteration = 0 ... max_nonlinear_iterations:
  export x_solution_ -> dofs_
  export dofs_ -> locally_relevant_dofs_
  run state updater callbacks
  assemble residual F(x)
  compute ||F||_2
  check convergence
  assemble sparse Jacobian or structured tree linearization
  solve J * delta = -F
  compute ||delta||_2
  x_solution_ = x_solution_ + delta
```

The convergence check requires both:

```text
||F||_2 <= nonlinear_residual_tolerance
||delta||_2 <= nonlinear_increment_tolerance
```

The increment check is treated as satisfied on iteration 0, matching the idea that there is no previous Newton increment before the first correction.

The critical integration point is:

```cpp
linear_solver_->linearization_type()
```

If the linear solver reports `SparseJacobian`, `NewtonSolver` calls sparse Jacobian assemblers and completes the sparse matrix if necessary.

If the linear solver reports `StructuredTreeBlocks`, `NewtonSolver` resets or clears `TreeLinearization`, calls all tree-linearization assemblers, and then passes the result into the tree linear solver with:

```cpp
linear_solver_->set_tree_linearization(tree_linearization_);
```

This is how the tree solver is connected to the custom Newton loop.

## Tree Metadata

The tree metadata is built by:

```cpp
build_reduced_lung_tree_metadata(ReducedLungTreeMetadataContext{...})
```

It is built from existing reduced-lung setup data, not from a new input format.

Inputs include:

- `ReducedLungParameters`, especially `lung_tree.topology`.
- `first_global_dof_of_ele` and `global_dof_per_ele`.
- Airway and terminal-unit model containers.
- Connection and bifurcation containers.
- Boundary-condition containers.
- The row map and locally relevant dof map.

### Topology Representation

The input topology gives each element two nodes:

```text
element_nodes[0] = inlet / parent side
element_nodes[1] = outlet / child side
```

The metadata builder uses this directed connectivity to find parent-child relations. It builds a map of elements starting at each node. A parent element's children are the elements whose inlet node equals the parent outlet node.

Each `TreeElementMetadata` stores:

- `global_element_id`.
- `kind`, either `Airway` or `TerminalUnit`.
- `inlet_node_id` and `outlet_node_id`.
- `parent_element_index`.
- up to two `child_element_indices`.
- `child_count`.
- first global dof and number of dofs.
- global dof ids and local dof ids.
- first local and global state-equation ids.
- number of state equations.
- owner MPI rank.

The metadata also stores:

- one `root_element_index`.
- one `root_node_id`.
- airway and terminal-unit element-index lists.
- `top_down_layers`.
- `bottom_up_layers`, which are the top-down layers reversed.
- global dof count, global equation count, and locally relevant dof count.

### Junction Metadata

For internal tree nodes, the reduced-lung formulation creates either a connection or a bifurcation.

A connection has one child and two equations:

```text
p_out_parent - p_in_child = 0
q_out_parent - q_in_child = 0
```

A bifurcation has two children and three equations:

```text
p_out_parent - p_in_child_1 = 0
p_out_parent - p_in_child_2 = 0
q_out_parent - q_in_child_1 - q_in_child_2 = 0
```

`TreeJunctionMetadata` stores the junction kind, parent and child element indices, equation ids, owner rank, and the involved dofs.

### Boundary Metadata

Boundary rows are stored in `TreeBoundaryConditionMetadata`.

Each boundary metadata entry records:

- boundary type, pressure or flow.
- side, inlet or outlet.
- attached node id.
- attached element index.
- local and global equation id.
- constrained global and local dof id.
- owner MPI rank.

The tree solver requires a root inlet boundary and one outlet boundary for every leaf element. In the current root closure implementation, the root inlet boundary effectively needs to constrain the root inlet pressure dof. A root flow closure would require extra logic that combines the root subtree `(G, h)` relation with the root boundary row.

### Metadata Validation

The metadata builder checks the structural assumptions required by the tree algorithm:

- At least one element exists.
- Every element has exactly two valid topology nodes.
- Element endpoints are not identical.
- Each element has known dof offset and dof count.
- Each element has model equation metadata.
- Each element has at least three dofs.
- The global equation count equals the global dof count.
- Directed children have at most one parent.
- Branch degree is at most two.
- The topology is acyclic.
- Exactly one root exists.
- All elements are reachable from the root.
- Terminal units are leaves.
- Junction metadata matches the directed topology.
- Every one-child parent has connection metadata.
- Every two-child parent has bifurcation metadata.
- Every boundary condition attaches to the expected element side.
- Every boundary condition constrains the expected dof for its side and type.
- The root inlet side has a boundary condition.
- Every leaf outlet side has a boundary condition.

These checks are safety and validation checks, not performance optimizations.

## Unknowns, Rows, Dofs, Elements, Junctions, and Boundary Conditions

The tree solver works on the current duplicated-endpoint dof layout.

For each element, the element-local unknown vector used in bottom-up elimination is:

```text
u = all element dof corrections except inlet pressure correction
```

The inlet pressure correction is not eliminated locally. It is the scalar parameter for the subtree relation.

The usual element dof layouts are:

| Element type | Dofs | State equations | Tree block size |
| --- | --- | ---: | ---: |
| Rigid airway | `p1, p2, q1` | 1 | 2 |
| Kelvin-Voigt airway | `p1, p2, q1, q2` | 2 | 3 |
| Terminal unit | `p1, p2, q` | 1 | 2 |

In all cases:

- `p1` is the inlet pressure dof.
- `p2` is the outlet pressure dof.
- `global_dof_ids[2]` is treated as the inlet-flow dof.
- The last element dof is treated as the outlet-flow dof in junction metadata.

For a leaf element, the local dense block rows are:

```text
element state-equation rows
leaf outlet boundary row
```

For an internal element, the local dense block rows are:

```text
element state-equation rows
downstream junction flow-conservation row
```

Pressure-continuity junction rows are not included directly in the dense local block. They are used separately to express child inlet pressure as an affine function of parent outlet pressure.

## Linearization / Coefficient Source

There are two coefficient-source modes in the serial tree solver:

```cpp
TreeNewtonLinearSolverCoefficientSource::SparseJacobian
TreeNewtonLinearSolverCoefficientSource::StructuredTreeBlocks
```

### Sparse Jacobian Source

The sparse source reads coefficients from a completed sparse Jacobian using `SparseMatrix::extract_my_row_view(...)`.

This source is still available and is important for validation tests because it answers this question:

```text
If the tree solver reads the exact sparse Jacobian entries, does it produce the same delta as UMFPACK?
```

The sparse source is not the current runtime `NewtonTree` path.

### Structured Tree Blocks Source

Runtime `NewtonTree` uses structured tree blocks.

The structured path is:

```text
assemble residual F
assemble TreeLinearization
TreeNewtonLinearSolver reads structured coefficients
bottom-up/top-down tree solve
write delta
```

`TreeLinearization` stores coefficients by local residual row and local locally-relevant dof id:

```text
rows_[local_row_id] = vector of (local_dof_id, value)
```

The structured tree-linearization callbacks are registered in `create_default_reduced_lung_assembly_pipeline(...)` in this order:

```text
Airways
Terminal units
Junctions
Boundary conditions
```

The derivative formulas remain in the physics modules. The tree solver consumes the resulting coefficients. This avoids duplicating physics formulas inside the solver.

### Current Coefficient Implementation Details

The serial structured solver now precomputes the coefficient locations it needs:

- root boundary coefficient.
- inlet-pressure coefficient for each local equation row.
- dense matrix coefficients for each local block.
- parent-pressure and child-pressure coefficients for pressure-continuity rows.
- child-flow coefficients for junction flow rows.

When a new `TreeLinearization` is attached, the solver resolves the row-entry indices and copies numeric values into direct coefficient buffers. During the serial structured solve, hot loops read these direct buffers instead of repeatedly searching `TreeLinearization` rows.

The distributed solver currently consumes structured tree blocks only. It maps global row and dof ids to local row and locally relevant dof ids on each rank, then reads the coefficient from the rank-local `TreeLinearization`.

### Current Runtime Summary

The current runtime behavior is:

| Workflow | Coefficient source |
| --- | --- |
| `Nox` | Sparse Jacobian through NOX. |
| `NewtonSparse` | Sparse Jacobian through `NewtonSolver`. |
| `NewtonTree`, serial | Structured `TreeLinearization`, no sparse Jacobian assembly. |
| `NewtonTree`, MPI | Structured `TreeLinearization`, no sparse Jacobian assembly. |
| Serial tree validation tests | Both sparse-Jacobian tree source and structured tree source. |

## Bottom-Up Algorithm

The serial algorithm is implemented in `TreeNewtonLinearSolver::solve(...)`. The distributed solver implements the same mathematical operations but only for rank-owned elements and exchanges boundary information between ranks.

### Local Element Relation

For one element, define:

```text
p = delta_p_in
u = corrections of all other element dofs
```

After child subtrees have been substituted, the element local system has the form:

```text
A * u = b_const + b_pin * p
```

where:

- `A` contains coefficients with respect to local unknowns `u`.
- `b_const` contains `-F` and child-subtree intercept contributions.
- `b_pin` is the negative coefficient of inlet pressure in each local row.

The solver solves two right-hand sides with the same dense local matrix:

```text
A * intercept = b_const
A * slope     = b_pin
```

Therefore every local unknown correction is represented as:

```text
u = intercept + slope * delta_p_in
```

The element subtree relation is the inlet-flow entry of this affine solution:

```text
G = slope[inlet_flow_unknown_index]
h = intercept[inlet_flow_unknown_index]
delta_q_in = G * delta_p_in + h
```

This `(G, h)` pair is stored for the element and used by its parent.

### Child Subtree Substitution

Assume a child subtree has already been condensed to:

```text
delta_q_child_in = G_child * delta_p_child_in + h_child
```

The pressure-continuity row between the parent and that child is treated as:

```text
a_parent * delta_p_parent_out + a_child * delta_p_child_in = rhs_pressure
```

where:

```text
rhs_pressure = -F_pressure_row
```

Solving for the child inlet pressure gives:

```text
delta_p_child_in = pressure_slope * delta_p_parent_out + pressure_intercept
pressure_slope = -a_parent / a_child
pressure_intercept = rhs_pressure / a_child
```

Then the child flow relation becomes:

```text
delta_q_child_in
  = G_child * (pressure_slope * delta_p_parent_out + pressure_intercept) + h_child
```

The parent flow-conservation row contains the child inlet flow coefficient. The solver adds the child-flow slope contribution into the parent local matrix entry for parent outlet pressure and subtracts the child-flow intercept contribution from the RHS.

This is the Schur-complement step. It removes the entire child subtree from the parent local solve and replaces it with one modified coefficient and one RHS shift.

### Dense Local Solve

The serial solver uses a small internal Gaussian-elimination routine with partial pivoting. It solves the intercept and slope right-hand sides together. If a pivot is below `pivot_tolerance`, the solver throws a singular or underconstrained local block error.

The default pivot tolerance is:

```text
1.0e-12
```

The current optimized serial path also has specialized `2x2` and `3x3` batch handling for common element block sizes. The mathematical result is the same as the generic dense solve.

### Bottom-Up Traversal

The solver processes:

```text
tree_metadata_.bottom_up_layers       # scalar path / small trees
bottom_up_layer_groups_               # grouped path for larger trees
```

from leaves to root.

After a layer is complete, every element in that layer has a valid subtree relation `(G, h)`. Parent layers can then use those child relations. The grouped traversal is derived from the same bottom-up layers; it only clusters elements by block size and child count so the hot solve can use shape-specific loops.

Implementation trace in `TreeNewtonLinearSolver::solve(...)`:

- `build_symbolic_plan()` has already precomputed element offsets, local row lists, matrix slots, child-interface slots, and grouped traversals.
- `assemble_scalar_element(...)` or the grouped assembly path fills `workspace_matrix_`, `workspace_rhs_constant_`, and `workspace_rhs_inlet_pressure_`.
- `solve_dense_system(...)`, `solve_2x2_batch(...)`, or `solve_3x3_batch(...)` fills `workspace_intercept_` and `workspace_slope_`.
- `write_subtree_relation(...)` stores `subtree_relation_G_` and `subtree_relation_h_` from the inlet-flow unknown.
- Child pressure-continuity rows store `child_pressure_slope_` and `child_pressure_intercept_` for the later top-down pass.

## Root Closure

After bottom-up condensation, the root element has its subtree relation. For the current runtime pressure-boundary cases, the root inlet pressure correction is determined directly by the root inlet boundary row.

The boundary linearized row is:

```text
boundary_coeff * delta_p_root = -F_boundary
```

In the implementation, `rhs_value(residual, row)` returns `-residual[row]`, so the code computes:

```text
delta_p_root = rhs_value(root_boundary_row) / root_boundary_coeff
```

The root subtree relation is then used in top-down recovery to compute root inlet flow and the other root element unknowns.

Important limitation:

```text
The current root closure is effectively a root pressure closure.
```

The metadata supports boundary-condition types, but `TreeNewtonLinearSolver` closes the root using the root inlet pressure dof. A root flow boundary would need a different closure equation that combines the root boundary row with the root subtree `(G, h)` relation.

## Top-Down Recovery

The top-down phase starts with the known root inlet pressure correction.

For each element, the solver already has:

```text
u = intercept + slope * delta_p_in
```

The top-down recovery does:

```text
known delta_p_in for current element
write delta_p_in into delta
for each local unknown:
  value = slope[unknown] * delta_p_in + intercept[unknown]
  write value into delta
if element has children:
  recover parent outlet pressure
  use stored pressure-continuity affine relation
  stamp child inlet pressure for next layer
```

The pressure-continuity relation saved during bottom-up is:

```text
delta_p_child_in = child_pressure_slope * delta_p_parent_out + child_pressure_intercept
```

The solver processes:

```text
tree_metadata_.top_down_layers        # scalar path / small trees
top_down_layer_groups_                # grouped path for larger trees
```

from root to leaves.

In the serial solver, inlet-pressure stamp arrays guard against missing or duplicate inlet-pressure writes. This is a correctness check, not a performance optimization.

## Output Correction Vector

The output of the tree solver is the same object expected by the custom Newton loop:

```text
Core::LinAlg::Vector<double>& delta
```

It uses the same global dof numbering as the sparse solver correction vector.

### Serial Output

The serial solver requires all correction dofs to be locally available. It validates this once and caches local correction indices.

During top-down recovery, it writes directly into local vector storage:

```text
delta_values[local_dof_id] = value
```

This replaced repeated global replacement calls and is a performance optimization.

The solver also validates during symbolic setup that each global correction dof is written exactly once by top-down recovery.

### Distributed Output

The distributed solver computes element corrections on the rank that owns the element. The correction vector itself may be distributed differently, so the distributed solver uses a final correction scatter:

```text
global_dof_id, delta_value
```

If the current rank owns the correction-vector entry, it writes directly. Otherwise it queues a scalar message to the correction owner. These messages are exchanged after top-down recovery.

## Relation to the Python/Jupyter Prototype

The relevant prototype is:

```text
info_oc/05_airway_tree_solver_kv_soa_numpy_nonlinear.ipynb
```

### What Is Mathematically the Same

The main idea is the same:

- Preprocess the tree into parent-child relations.
- Build bottom-up and top-down traversal layers.
- Avoid recursion in the repeated solve.
- Condense each subtree into `dQ = G dPin + h`.
- Use a root pressure correction to start recovery.
- Recover pressure and flow corrections from root to leaves.
- Store solver quantities in array-like layouts rather than object-heavy traversal during the hot solve.

### What Is Different in the C++ Implementation

The C++ implementation is more general and closer to the real 4C equations:

- It uses the existing duplicated endpoint dof formulation.
- It uses explicit junction rows instead of implicit child coupling formulas.
- It uses explicit boundary-condition rows.
- It uses existing local and global row/dof maps.
- It handles different element block sizes.
- It supports rigid airways, Kelvin-Voigt airways, nonlinear airway resistance derivatives, terminal-unit rheology, and terminal-unit elasticity through existing physics modules.
- It has both serial and MPI-capable tree solvers.
- It writes one global `delta` vector that the custom Newton loop applies.
- It compares directly against the sparse solver and NOX workflows.

The prototype had separate airway and terminal-unit arrays such as `aw_G`, `aw_h`, `tu_G`, and `tu_h`. The C++ solver stores one generic `G` and `h` relation per tree element. This is necessary because in 4C both airways and terminal units are elements in the same dof/row system.

The prototype hard-coded simplified model formulas. The C++ implementation deliberately keeps derivative formulas in the physics modules and reads the assembled coefficients. This avoids having two separate versions of the airway, terminal-unit, junction, and boundary derivatives.

### Why the C++ Implementation Was Done This Way

The chosen implementation is conservative and validation-friendly:

- It preserves the current equations, row order, and dof order.
- It keeps NOX and sparse Newton as reference workflows.
- It allows sparse-Jacobian tree coefficient extraction for direct validation.
- It allows runtime `NewtonTree` to avoid sparse Jacobian assembly through structured `TreeLinearization`.
- It can be tested one Newton correction at a time against the sparse solver.
- It can be integrated into the existing timestep loop without rewriting the reduced-lung physics.

## Current Optimizations

The current serial tree solver is already substantially optimized compared with the first validation implementation.

Implemented optimizations include:

- Runtime `NewtonTree` uses structured `TreeLinearization`, so it avoids sparse Jacobian assembly and sparse matrix completion.
- The serial solver builds a symbolic plan once in its constructor and reuses it across Newton corrections.
- Element unknown offsets, equation offsets, matrix offsets, and child-interface offsets are precomputed.
- Workspace storage is flattened into contiguous vectors.
- Subtree relations are stored as separate `G` and `h` arrays.
- Elements are grouped by layer, block size, and child count.
- Common `2x2` bottom-up assembly shapes are specialized for leaf, one-child, and two-child cases.
- `2x2` dense solves are performed directly from workspace arrays with determinant formulas and dense fallback.
- `3x3` groups use a specialized batch-oriented path with fallback to the generic dense solver.
- Structured coefficient locations are precomputed.
- Structured coefficient numeric values are copied into direct buffers for hot solve loops.
- `TreeLinearization` storage is reused when dimensions are unchanged.
- Cached structured entry indices are reused when possible.
- Small trees use scalar traversal to avoid grouped staging overhead.
- Larger trees use grouped traversal.
- Small top-down groups use scalar recovery.
- Common `2x2` top-down recovery groups avoid temporary staging buffers.
- Serial top-down writes use cached local vector indices instead of repeated global replacement calls.
- Detailed timers are gated so no-profile serial solves avoid unnecessary timer calls.
- Optional profile structs record solve time, bottom-up time, top-down time, dense solve counts, lookup counts, and distributed communication counters.

The distributed solver has also been improved from an earlier replicated MPI strategy. It now performs rank-local bottom-up and top-down work and communicates only partition-boundary relations, child inlet pressures, and final correction scatter values.

## Safety and Validation Checks Versus Performance Optimizations

It is useful to separate checks that make the solver safe from changes that make it faster.

### Safety and Validation Checks

These are mainly correctness checks:

- Metadata rejects invalid topology, cycles, multiple roots, disconnected trees, unsupported branch degree, missing junctions, missing root boundaries, and missing leaf outlet boundaries.
- Metadata checks that terminal units are leaves.
- Metadata checks that the system is square.
- Boundary metadata checks that boundary rows constrain the expected dof for the side and type.
- The serial solver checks `comm_size == 1`.
- The sparse-source tree path checks that the sparse Jacobian is completed.
- The structured-source tree path checks that `TreeLinearization` exists and has expected dimensions.
- The solver asserts that residual rows and correction dofs are locally available where required.
- Required coefficients are checked for missing or near-zero values.
- Dense local solves check pivot magnitude against `pivot_tolerance`.
- Top-down recovery checks that inlet pressures are known before use.
- Stamp checks reject duplicate child inlet-pressure writes.
- Symbolic setup checks each correction dof is written exactly once in serial recovery.
- Distributed message receivers validate incoming element and dof ids.

### Performance Optimizations

These are mainly speed and memory-layout changes:

- Structured `TreeLinearization` avoids sparse Jacobian assembly for runtime `NewtonTree`.
- Symbolic plan reuse avoids rebuilding offsets and traversal plans.
- Flattened arrays improve locality.
- Same-shape layer groups allow batch-style loops.
- Direct coefficient buffers avoid row-entry lookup in the structured serial solve loop.
- Specialized `2x2` and `3x3` paths reduce dense-solve overhead.
- Direct local `delta` writes avoid repeated global-id lookup in serial top-down recovery.
- Timer gating avoids profiling overhead when no profile is attached.
- Distributed partition-boundary messaging avoids gathering all coefficients and residuals to every rank.

Some mechanisms serve both purposes. For example, dense fallback for near-singular `2x2` blocks preserves robustness while the fast determinant path improves performance.

## Remaining Optimization Opportunities

### Avoiding Sparse Jacobian Assembly

Runtime `NewtonTree` already avoids sparse Jacobian assembly. Sparse Jacobian support remains for `Nox`, `NewtonSparse`, and validation tests.

Future work should ensure all production tree workflows stay on the structured path and avoid accidentally requiring sparse matrix completion for diagnostics or fallback behavior.

### Structured Local Derivative Blocks

`TreeLinearization` is still a generic row-oriented local coefficient container. It is lighter than a sparse matrix, but it is not yet a fully preplanned dense block storage format.

A faster design would let physics modules write directly into precomputed element, junction, and boundary coefficient slots used by the tree solver.

Possible direction:

- Keep derivative formulas in physics modules.
- Replace generic `set_value(row, dof, value)` calls with structured block writer callbacks.
- Precompute exact storage slots during the tree solver symbolic setup.
- Write derivative values directly into tree-solver numeric arrays.
- Remove remaining row/dof searches from structured assembly.

### Reducing Coefficient Extraction Overhead

Serial structured solve-time coefficient lookup has mostly been removed from the hot loop. Remaining overhead is now shifted into structured assembly and coefficient refresh.

Potential follow-ups:

- Cache more writer locations in physics modules.
- Avoid repeated `TreeLinearization::set_value` row scans.
- Use fixed-size arrays for common rows.
- Separate symbolic sparsity from numeric values in `TreeLinearization`.

### Improving Memory Layout

The serial solver already uses flattened workspaces and grouped traversal. Further improvements could include:

- More direct dense-block storage by element shape.
- Specialized closed-form `3x3` kernels if benchmark evidence supports it.
- Less staging between assembly and solve arrays.
- Better cache-aware ordering for very large trees.
- Release-build threshold tuning for scalar versus grouped traversal.

### Improving Parallel and MPI Support

The distributed solver now communicates only partition-boundary data, but it still uses collective exchanges per tree layer.

Potential improvements:

- Cache correction-owner lookup instead of querying it every solve.
- Precompute per-layer send and receive plans from metadata.
- Replace `MPI_Alltoallv` layer exchanges with point-to-point communication between ranks that actually share tree cut edges.
- Overlap local condensation or recovery with communication.
- Add NP4 and larger distributed benchmarks.
- Add distributed mixed airway/terminal-unit validation cases.
- Add distributed full-workflow tests, not only correction tests.

### Improving Profiling and Diagnostics

Existing profiling is useful, but further diagnostics would help thesis reporting and debugging:

- Release-build benchmark runs with stable CPU settings.
- Allocation profiling for tree assembly and solve phases.
- Separate timing for residual assembly, structured assembly, coefficient refresh, bottom-up solve, top-down recovery, and MPI communication.
- Optional tree linear residual checks against `J * delta + F` for small validation cases.
- Pivot statistics and worst-local-block diagnostics.
- Per-layer element counts, group sizes, and communication sizes.

## Tests and Benchmarks That Validate the Tree Solver

The current validation strategy compares the tree solver against sparse and NOX paths at several levels.

### Serial Correction Tests

File:

```text
src/reduced_lung/tests/4C_reduced_lung_tree_linear_solver_test.cpp
```

These tests assemble a residual and sparse Jacobian, assemble structured tree linearization, and compare corrections from:

```text
SparseNewtonLinearSolver
sparse-source TreeNewtonLinearSolver
structured-source TreeNewtonLinearSolver
```

The correction vector is compared entry-by-entry.

Covered cases include:

- Single terminal unit.
- Serial rigid airways.
- Rigid-airway bifurcation.
- Kelvin-Voigt airways.
- Nonlinear airway resistance.
- Four-element Maxwell terminal unit.
- Mixed airway and terminal-unit tree.
- Reusing the same structured tree solver across several states and times.

### Full Workflow Tests

The same test file compares complete time-step workflows for:

```text
Nox
NewtonSparse
NewtonTree
```

Covered cases include:

- Single terminal unit.
- Serial rigid airways.
- Bifurcation rigid airways.
- Mixed airways and terminal units.

The tests compare solution vectors, owned dofs, residual norms, flow balances, and terminal-unit volume updates where relevant.

### Distributed Correction Tests

File:

```text
src/reduced_lung/tests/4C_reduced_lung_distributed_tree_linear_solver_test.np2.cpp
```

These two-rank tests compare:

```text
SparseNewtonLinearSolver correction
DistributedTreeNewtonLinearSolver correction
```

Covered cases are:

- Serial airway chain on two ranks.
- Bifurcation airways on two ranks.

Both tests assert that the metadata contains at least one cross-rank edge, so partition-boundary communication is actually exercised.

### Metadata Tests

File:

```text
src/reduced_lung/tests/4C_reduced_lung_tree_metadata_test.cpp
```

These tests validate tree metadata construction, including connection metadata, bifurcation metadata, traversal layers, terminal-unit leaves, cycle rejection, unsupported branch-degree rejection, and missing root boundary rejection.

### Runtime Input Tests

The implementation notes record runtime YAML variants for explicit `NewtonSparse` and `NewtonTree` selection. These check that the executable path selects the intended workflow.

The large generation-10 airways-only runtime inputs in `Implementation_NextStep_9.md` are meant for solver-path timing and stress testing, not physiological validation.

### Benchmarks

Benchmark files:

```text
src/reduced_lung/benchmark_tests/4C_reduced_lung_solver_benchmark.cpp
src/reduced_lung/benchmark_tests/4C_reduced_lung_distributed_tree_solver_benchmark.cpp
```

Serial benchmark groups include:

- Full solve comparisons for `Nox`, `NewtonSparse`, and `NewtonTree`.
- Assembly-phase timing.
- Sparse linear solve timing.
- Structured tree linear solve timing.

Distributed benchmark groups include:

- Distributed full solve comparisons.
- Distributed structured tree linear solve timing.

The historical optimization notes record release benchmark measurements, but this document did not rerun tests or benchmarks.
