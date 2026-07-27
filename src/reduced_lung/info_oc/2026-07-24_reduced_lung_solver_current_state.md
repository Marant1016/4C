# Reduced-Lung Solver Current State - 2026-07-24

This document summarizes the current non-1D `src/reduced_lung/` solver implementation after the custom Newton, structured tree, distributed tree, validation, and profiling work. It is written for a reader who knows the original reduced-lung code and its NOX-based workflow, but has not followed the new implementation phases.

The `src/reduced_lung/src/1d_pipe_flow/` path is separate and is not described here.

## Executive Summary

The reduced-lung solver now has three selectable nonlinear solver workflows:

```text
Nox
NewtonSparse
NewtonTree
```

`Nox` remains the default. Existing YAML files that omit `reduced_dimensional_lung.dynamics.nonlinear_solver` keep the original NOX behavior.

`NewtonSparse` and `NewtonTree` are explicit opt-in workflows. Both use the new custom full-step Newton driver. `NewtonSparse` solves each Newton correction with the existing sparse `Core::LinAlg::Solver` backend. `NewtonTree` solves each Newton correction with a tree-structured linear solver. In the current runtime implementation, `NewtonTree` assembles structured derivative blocks into `TreeLinearization` and does not assemble or complete the global sparse Jacobian for Newton corrections.

The original physics formulation is preserved. Element residuals, junction equations, boundary-condition equations, dof ordering, and row ordering are still the reduced-lung duplicated-endpoint formulation. The tree solver exploits that existing structure instead of changing the equations.

## Original NOX-Based Workflow

The original reduced-lung workflow solved one nonlinear algebraic system per time step:

```text
F(x, time) = 0
```

The unknown vector `x` contains element-local endpoint pressures and flows. Elements do not share nodal unknowns. Instead, each element owns its endpoint dofs and junction equations enforce pressure continuity and flow conservation between neighboring elements.

The original runtime path was:

```text
ReducedLung::reduced_lung_main()
  -> ReducedLungSimulation::initialize()
  -> ReducedLungSimulation::run()
  -> ReducedLungSimulation::solve_timestep()
  -> NoxSolver::solve(current_time)
  -> NOX::Nln::Adapter::solve()
```

In that workflow, reduced lung supplied NOX with:

- the current nonlinear solution vector,
- a residual callback,
- a sparse Jacobian callback,
- a sparse matrix operator,
- the selected 4C linear solver from `dynamics.linear_solver`, normally UMFPACK in the regression inputs,
- absolute residual and increment tolerances,
- the maximum nonlinear iteration count.

NOX owned the nonlinear iteration loop. It called the reduced-lung residual and Jacobian callbacks, solved the Newton linear system through 4C linear-solver infrastructure, checked convergence, applied the update, and copied the final solution back into the bound reduced-lung vector.

The NOX parameter list still uses:

- `Nonlinear Solver = Line Search Based`,
- Newton direction,
- `Line Search/Method = Full Step`,
- absolute two-norm residual convergence,
- absolute update convergence with first update skipped,
- maximum nonlinear iteration guard.

Important sign detail in the NOX path: the NOX linear-system layer solves an internal `J * y = F` system and the NOX Newton direction machinery applies the negative sign afterward. Mathematically, the resulting nonlinear correction is equivalent to:

```text
J * delta = -F
x_new = x_old + delta
```

## Important Architectural Changes

The new implementation keeps the NOX path and adds custom solver infrastructure beside it.

### Assembly Pipeline Generalization

The old NOX-specific callback grouping was generalized to:

```text
ReducedLungAssemblyPipeline
```

Implemented in:

```text
src/reduced_lung/src/4C_reduced_lung_helpers.hpp
src/reduced_lung/src/4C_reduced_lung_helpers.cpp
```

The pipeline contains four callback groups:

- residual assemblers,
- sparse Jacobian assemblers,
- structured tree-linearization assemblers,
- nonlinear state updaters.

The compatibility alias remains:

```text
using NoxAssemblyPipeline = ReducedLungAssemblyPipeline;
```

The callback order is the same for residuals, sparse Jacobians, and structured tree linearizations:

```text
Airways
Terminal units
Junctions
Boundary conditions
```

The state updater order is:

```text
Airways
Terminal units
```

### Custom Newton Solver

The custom nonlinear solver is:

```text
NewtonSolver
NewtonSolverContext
```

Implemented in:

```text
src/reduced_lung/src/4C_reduced_lung_newton_solver.hpp
src/reduced_lung/src/4C_reduced_lung_newton_solver.cpp
```

`NewtonSolver` owns a full-step Newton iteration loop. It reuses the same residual and state-update callbacks as NOX. Depending on the selected linear solver, it assembles either a sparse Jacobian or structured tree-linearization blocks.

### Newton Linear Solver Interface

The Newton correction solve is abstracted behind:

```text
NewtonLinearSolver
NewtonLinearSystemMetadata
NewtonLinearizationType
```

Implemented in:

```text
src/reduced_lung/src/4C_reduced_lung_linear_solver.hpp
src/reduced_lung/src/4C_reduced_lung_linear_solver.cpp
```

The interface operation is:

```text
solve(jacobian, residual, x, metadata, delta)
```

The documented convention is:

```text
jacobian * delta = -residual
x_new = x_old + delta
```

The interface can report whether it needs:

```text
SparseJacobian
StructuredTreeBlocks
```

This lets `NewtonSolver` assemble the correct linearization source for each correction.

### Sparse Newton Linear Solver

The sparse backend is:

```text
SparseNewtonLinearSolver
SparseNewtonLinearSolverContext
```

It wraps `Core::LinAlg::Solver`, forms `rhs = -residual`, zeros `delta`, and calls:

```text
Core::LinAlg::Solver::solve(jacobian, delta, rhs, solver_params)
```

It sets:

```text
refactor = true
reset = metadata.nonlinear_iteration == 0
```

This backend is used by the runtime `NewtonSparse` workflow and by many tests as the reference for tree corrections.

### Tree Metadata

Explicit tree/layout metadata is built by:

```text
build_reduced_lung_tree_metadata(context)
```

Implemented in:

```text
src/reduced_lung/src/4C_reduced_lung_tree_metadata.hpp
src/reduced_lung/src/4C_reduced_lung_tree_metadata.cpp
```

The tree metadata is built from existing reduced-lung setup data. It does not introduce a new input format.

### Structured Tree Linearization

The current runtime `NewtonTree` path uses:

```text
TreeLinearization
```

Implemented in:

```text
src/reduced_lung/src/4C_reduced_lung_tree_linearization.hpp
src/reduced_lung/src/4C_reduced_lung_tree_linearization.cpp
```

`TreeLinearization` stores derivative coefficients addressed by:

```text
local residual row id
local locally-relevant dof id
```

This mirrors the coefficient lookup that the early validation tree solver performed on a completed sparse Jacobian, but it avoids constructing the generic sparse matrix for runtime `NewtonTree` corrections.

### Serial And Distributed Tree Linear Solvers

The serial tree solver is:

```text
TreeNewtonLinearSolver
```

The distributed tree solver is:

```text
DistributedTreeNewtonLinearSolver
```

Declared in:

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp
```

Implemented in:

```text
src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp
src/reduced_lung/src/4C_reduced_lung_distributed_tree_linear_solver.cpp
```

The serial solver can consume either sparse Jacobian coefficients or structured tree blocks. Runtime `NewtonTree` constructs it with `StructuredTreeBlocks` when `comm_size == 1`.

The distributed solver consumes structured tree blocks only. Runtime `NewtonTree` constructs it when `comm_size > 1`.

### Profiling And Benchmarking

Optional profiling structs were added in:

```text
src/reduced_lung/src/4C_reduced_lung_solver_profile.hpp
```

A reduced-lung Google Benchmark harness was added in:

```text
src/reduced_lung/benchmark_tests/4C_reduced_lung_solver_benchmark.cpp
src/reduced_lung/benchmark_tests/CMakeLists.txt
```

The profile pointers default to `nullptr`, so normal runtime behavior is unchanged unless benchmark code attaches profiles.

## Current Runtime Solver Options

The solver workflow is selected through:

```yaml
reduced_dimensional_lung:
  dynamics:
    nonlinear_solver: Nox
```

The C++ enum is:

```text
ReducedLungParameters::NonlinearSolverType
```

Defined in:

```text
src/reduced_lung/src/4C_reduced_lung_input.hpp
```

Registered in valid input parameters in:

```text
src/reduced_lung/src/4C_reduced_lung_input.cpp
```

The supported values are:

| Option | Nonlinear driver | Linear correction solver | Current role |
| --- | --- | --- | --- |
| `Nox` | NOX through `NoxSolver` | `Core::LinAlg::Solver` through NOX | Default/reference workflow |
| `NewtonSparse` | `NewtonSolver` | `SparseNewtonLinearSolver` | Custom Newton sparse reference/debug workflow |
| `NewtonTree` | `NewtonSolver` | `TreeNewtonLinearSolver` or `DistributedTreeNewtonLinearSolver` | Structured tree-solver workflow |

The default is intentionally:

```text
Nox
```

`NewtonSparse` and `NewtonTree` are explicit opt-in values. Existing YAML decks that omit `nonlinear_solver` preserve the baseline NOX behavior.

## Complete Current Workflow

### Input And Context Construction

Runtime starts in:

```text
ReducedLung::reduced_lung_main(Global::Problem& problem)
```

Implemented in:

```text
src/reduced_lung/src/4C_reduced_lung_main.cpp
```

`make_reduced_lung_context_from_problem(problem)` collects:

- `ReducedLungParameters` from `reduced_dimensional_lung`,
- the local MPI communicator,
- rebalance parameters,
- IO parameters,
- linear solver parameters selected by `dynamics.linear_solver`,
- the solver parameter callback,
- output control,
- the function manager used by time-dependent boundary conditions.

Then `ReducedLungSimulation` is constructed, initialized, and run.

### Initialization

`ReducedLungSimulation::initialize()` performs:

```text
validate_parameters()
build_discretization()
build_element_models()
build_node_entities()
assign_equation_ids()
build_maps_and_local_ids()
build_linear_system_and_solver()
```

Validation checks positive `time_increment`, non-negative `number_of_steps`, positive `results_every`, and positive `max_nonlinear_iterations`.

`build_discretization()` builds a lightweight line-element `Core::FE::Discretization` from the input topology. Input node and element IDs are treated as 1-based in the YAML fields and converted to zero-based internal IDs. Element node order is directional:

```text
element_nodes[0] = inlet / parent side
element_nodes[1] = outlet / child side
```

`build_element_models()` creates local airway and terminal-unit model blocks, computes global dof offsets, assigns global dof IDs, and attaches evaluators.

`build_node_entities()` creates the global node-to-element adjacency map, boundary-condition models, junction data, and evaluator callbacks.

`assign_equation_ids()` assigns local equation rows in this order:

```text
airway state equations
terminal-unit state equations
junction equations
boundary-condition equations
```

`build_maps_and_local_ids()` creates:

- locally owned dof map,
- row/equation map,
- locally relevant dof/column map.

It then assigns global equation IDs and local dof IDs to junction and boundary data, and local dof IDs to airway and terminal-unit data.

### Solver Construction

`build_linear_system_and_solver()` creates:

- `dofs_` on the locally owned dof map,
- `locally_relevant_dofs_` on the column map,
- `x_` on the row map,
- `sysmat_` as a sparse matrix with row map and locally relevant dof map,
- `ReducedLungAssemblyPipeline`.

The sparse matrix is allocated with an estimate of 4 entries per row. This supports Kelvin-Voigt airway rows that insert 4 coefficients.

Then it switches on `dynamics.nonlinear_solver`:

```text
Nox
  -> build_nox_solver()

NewtonSparse
  -> build_newton_solver_with_sparse_linear_solver()
  -> build_newton_solver()

NewtonTree
  -> build_newton_solver_with_tree_linear_solver()
  -> build_reduced_lung_tree_metadata(...)
  -> comm_size == 1 ? TreeNewtonLinearSolver : DistributedTreeNewtonLinearSolver
  -> build_newton_solver()
```

Even for `NewtonTree`, the common `sysmat_` object still exists because it is part of the linear-solver interface. Runtime `NewtonTree` does not assemble the sparse Jacobian because the selected tree linear solver reports `StructuredTreeBlocks` as its required linearization type.

### Time Integration

`ReducedLungSimulation::run()` loops over time steps. Each step calls:

```text
solve_timestep(step)
write_output_if_due(step)
```

`solve_timestep(step)` does:

```text
current_time_ += dt_
selected nonlinear solve at current_time_
TerminalUnits::end_of_timestep_routine(...)
Airways::end_of_timestep_routine(...)
```

The selected nonlinear solve is either:

```text
nox_solver_->solve(current_time_)
```

or:

```text
newton_solver_->solve(current_time_)
```

End-of-timestep updates are unchanged from the baseline. Terminal units update volume and rheology history. Airways update previous-step pressure/flow and wall-model history.

Runtime output is collected from `locally_relevant_dofs_` and written through `Core::IO::DiscretizationVisualizationWriterMesh` when `step % results_every == 0`.

## Custom Newton Algorithm

The custom Newton algorithm is implemented in:

```text
NewtonSolver::solve(double time)
```

The mathematical convention is always:

```text
J * delta = -F
x_new = x_old + delta
```

### Nonlinear Loop

For each time step, `NewtonSolver` performs iterations from 0 through `max_nonlinear_iterations`.

Each nonlinear iteration does:

```text
sync_state_from_x(x_solution_)
assemble residual F(x)
compute residual norm
check residual and increment convergence
assemble required linearization
solve for delta
compute increment norm
x_solution_ = x_solution_ + delta
```

If the initial residual is already converged, the solver can return with zero Newton corrections. The increment check is skipped on iteration 0 by treating it as converged, matching the intent of NOX's update-skip-first-iteration behavior.

### State Synchronization

`sync_state_from_x(x)` does:

```text
Core::LinAlg::export_to(x, dofs_)
Core::LinAlg::export_to(dofs_, locally_relevant_dofs_)
run state updater callbacks
```

The state updaters refresh model-internal nonlinear state before residual or derivative assembly. For example, nonlinear airway resistance and Kelvin-Voigt wall quantities can depend on the current trial dofs.

### Residual Assembly

`assemble_residual_for_current_state()` zeros the residual vector and runs the residual callbacks:

```text
Airways::update_residual_vector(...)
TerminalUnits::update_residual_vector(...)
Junctions::update_residual_vector(...)
BoundaryConditions::update_residual_vector(...)
```

The residual entries are the same equations as in the NOX path:

- airway pressure/flow and wall equations,
- terminal-unit rheology/elasticity equations,
- connection or bifurcation pressure-continuity equations,
- connection or bifurcation flow-conservation equations,
- boundary equations of the form `dof_value - prescribed_value`.

The residual norm is computed with the vector two-norm.

### Linearization Assembly

The linearization assembled by `NewtonSolver` depends on the selected `NewtonLinearSolver`.

For `NewtonSparse`:

```text
linear_solver_->linearization_type() == SparseJacobian
```

`NewtonSolver` runs the sparse Jacobian callbacks:

```text
Airways::update_jacobian(...)
TerminalUnits::update_jacobian(...)
Junctions::update_jacobian(...)
BoundaryConditions::update_jacobian(...)
```

It completes the sparse matrix if it is not already filled. The existing insert/replace sparse-matrix convention is preserved: on the first assembly, model blocks insert entries; on later assemblies, nonlinear coefficient values are replaced where needed, while constant junction and boundary blocks skip work once the matrix is filled.

For runtime `NewtonTree`:

```text
linear_solver_->linearization_type() == StructuredTreeBlocks
```

`NewtonSolver` resets `TreeLinearization` and runs:

```text
Airways::update_tree_linearization(...)
TerminalUnits::update_tree_linearization(...)
Junctions::update_tree_linearization(...)
BoundaryConditions::update_tree_linearization(...)
```

It then attaches the resulting `TreeLinearization` to the tree linear solver with:

```text
linear_solver_->set_tree_linearization(tree_linearization_)
```

No sparse Jacobian assembly or sparse matrix completion is performed by runtime `NewtonTree`.

### Linear Correction Solve

`solve_linear_correction(iteration)` builds:

```text
NewtonLinearSystemMetadata{
  current_time,
  time_step_size_dt,
  nonlinear_iteration
}
```

Then it calls the selected `NewtonLinearSolver`.

For `NewtonSparse`, `SparseNewtonLinearSolver` forms:

```text
rhs = -residual
```

and solves:

```text
J * delta = rhs
```

For `NewtonTree`, the serial or distributed tree solver uses the residual through the same convention. Internally, helper functions read each local RHS as:

```text
rhs_row = -residual_row
```

### Nonlinear Update

After the linear solver returns, `NewtonSolver` computes:

```text
increment_norm = ||delta||_2
```

Then it applies the full Newton step:

```text
x_solution_.update(1.0, delta_, 1.0)
```

This is:

```text
x_new = x_old + delta
```

There is currently no damping or line search in `NewtonSolver`.

### Convergence And Failure

The custom Newton convergence check requires both:

```text
||F||_2 <= nonlinear_residual_tolerance
||delta||_2 <= nonlinear_increment_tolerance
```

The increment condition is treated as satisfied on iteration 0. If convergence has not occurred by `max_nonlinear_iterations`, `NewtonSolver` throws with the current time, maximum correction count, final residual norm, and final increment norm.

## Tree-Based Solver: Mathematical Algorithm

The tree linear solver solves one Newton correction system:

```text
J * delta = -F
```

It exploits the directed tree topology rather than factoring the global sparse matrix.

The current formulation still uses duplicated endpoint dofs. Junction rows connect those duplicated dofs. The tree solver uses that formulation directly.

For each element subtree, the bottom-up pass condenses all downstream unknowns into a scalar affine relation at the element inlet:

```text
delta_q_in = G * delta_p_in + h
```

Here:

- `delta_p_in` is the inlet pressure correction for the root of the subtree,
- `delta_q_in` is the inlet flow correction of the same element,
- `G` is the condensed subtree conductance-like slope,
- `h` is the residual-dependent intercept.

This is a local Schur-complement elimination. Each child subtree is replaced by a scalar relation, then the parent eliminates its own local unknowns except inlet pressure. Repeating this from leaves to root condenses the whole tree.

Once the root subtree relation is available, the root inlet boundary condition determines the root inlet pressure correction. A top-down pass then recovers all element-local corrections and writes them into the global Newton correction vector.

## Tree Metadata

`ReducedLungTreeMetadata` is the structural input to the tree solvers. It records topology, ownership, dof layout, row layout, traversal layers, and boundary/junction closure information.

For each element, `TreeElementMetadata` stores:

- global element id,
- element kind, `Airway` or `TerminalUnit`,
- inlet and outlet node ids,
- parent element index,
- up to two child element indices,
- child count,
- first global dof,
- number of dofs,
- global dof ids,
- local dof ids in the locally relevant dof map,
- first local and global state-equation id,
- number of state equations,
- owner MPI rank.

For each junction, `TreeJunctionMetadata` stores:

- kind, `Connection` or `Bifurcation`,
- parent element index,
- child element indices,
- child count,
- first local and global equation id,
- number of junction equations,
- owner MPI rank,
- global and local dof ids involved in the junction block.

For each boundary condition, `TreeBoundaryConditionMetadata` stores:

- boundary type, pressure or flow,
- boundary side, inlet or outlet,
- node id,
- element index,
- local and global equation id,
- global and local constrained dof id,
- owner MPI rank.

The metadata also stores:

- root element index,
- root node id,
- airway element indices,
- terminal-unit element indices,
- bottom-up traversal layers,
- top-down traversal layers,
- global dof count,
- global equation count,
- locally relevant dof count.

The builder validates the assumptions required by the tree solver:

- topology has at least one element,
- every element has exactly two nodes,
- node ids are valid 1-based input ids,
- no self-connected element exists,
- every element has known dof offset and dof count,
- every element has model equation metadata,
- every element has at least three dofs,
- global equation count equals global dof count,
- every directed child has at most one parent,
- branch degree is at most two children,
- directed topology is acyclic,
- exactly one root element exists,
- all elements are reachable from the root,
- terminal units are leaves,
- connection and bifurcation metadata match the directed topology,
- every parent with one child has a connection,
- every parent with two children has a bifurcation,
- every boundary condition is attached to the inlet or outlet side of its element,
- boundary conditions constrain the expected dof for their side and type,
- the root inlet side has a boundary condition,
- every leaf outlet side has a boundary condition.

In MPI runs, metadata construction all-reduces local element equation metadata, junction metadata, and boundary metadata so every rank can build the same global directed tree. Local row and dof IDs remain valid only on ranks where those rows/dofs are locally available; they may be `-1` elsewhere.

## Element, Junction, And Boundary Blocks

### Element Blocks

Element blocks come from the existing airway and terminal-unit physics modules.

Airway files:

```text
src/reduced_lung/src/airways/4C_reduced_lung_airways.cpp
src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.cpp
src/reduced_lung/src/airways/4C_reduced_lung_airways_flow_resistance.cpp
```

Terminal-unit files:

```text
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit.cpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp
src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_elasticity.cpp
```

Rigid-wall airways have three dofs and one state equation:

```text
p1, p2, q1
```

Kelvin-Voigt airways have four dofs and two state equations:

```text
p1, p2, q1, q2
```

Terminal units have three dofs and one state equation:

```text
p1, p2, q
```

The structured tree-linearization callbacks write the same derivative coefficients as the sparse Jacobian callbacks, but into `TreeLinearization` with local row and local dof IDs.

Examples of structured element coefficients:

- rigid airway row: pressure coefficients on `p1` and `p2`, flow derivative on `q1`,
- Kelvin-Voigt airway momentum row: coefficients on `p1`, `p2`, `q1`, `q2`,
- Kelvin-Voigt airway mass/wall row: coefficients on `p1`, `p2`, `q1`, `q2`,
- terminal-unit row: coefficients on `p1`, `p2`, `q`, including elasticity and rheology derivatives.

### Junction Blocks

Junction blocks are implemented in:

```text
src/reduced_lung/src/4C_reduced_lung_junctions.hpp
src/reduced_lung/src/4C_reduced_lung_junctions.cpp
```

A connection has two equations:

```text
p_out_parent - p_in_child = 0
q_out_parent - q_in_child = 0
```

A bifurcation has three equations:

```text
p_out_parent - p_in_child_1 = 0
p_out_parent - p_in_child_2 = 0
q_out_parent - q_in_child_1 - q_in_child_2 = 0
```

The structured junction linearization writes constant coefficients into `TreeLinearization`. Pressure-continuity rows are later used during tree recovery to express child inlet pressures in terms of parent outlet pressure. Flow-conservation rows are used in the parent element's condensed local system.

### Boundary-Condition Blocks

Boundary-condition blocks are implemented in:

```text
src/reduced_lung/src/4C_reduced_lung_boundary_conditions.hpp
src/reduced_lung/src/4C_reduced_lung_boundary_conditions.cpp
```

Boundary residuals have the form:

```text
constrained_dof_value - prescribed_value = 0
```

The sparse and structured linearizations both insert a diagonal coefficient:

```text
dF / d(constrained_dof) = 1
```

The root inlet boundary row closes the root pressure correction. Leaf outlet boundary rows close leaf local systems during bottom-up condensation.

## Structured `TreeLinearization`

`TreeLinearization` is a row-oriented local coefficient container:

```text
rows_[local_row_id] = vector of (local_dof_id, value)
```

Its main operations are:

```text
reset(num_rows, num_dofs)
set_value(local_row_id, local_dof_id, value)
value(local_row_id, local_dof_id)
entries(local_row_id)
```

`set_value` updates an existing entry if one exists for the same local row and dof, otherwise it appends a new pair. `value` returns 0.0 if a coefficient is absent.

Runtime `NewtonTree` assembles one `TreeLinearization` per Newton correction. The linearization dimensions are:

```text
num_rows = local residual row count
num_dofs = local locally-relevant dof count
```

The current runtime path is therefore:

```text
sync nonlinear state
assemble residual F
assemble structured TreeLinearization
tree linear solve
update x
```

Sparse coefficient extraction still exists in `TreeNewtonLinearSolver` as `SparseJacobian` coefficient source. It is retained for serial sparse-versus-tree comparison tests, not for runtime `NewtonTree` construction.

## Serial Tree Solver: C++ Implementation And Data Flow

The serial solver class is:

```text
TreeNewtonLinearSolver
```

It is constructed with:

```text
TreeNewtonLinearSolverContext{
  tree_metadata,
  pivot_tolerance,
  coefficient_source,
  profile
}
```

For runtime `NewtonTree` in serial, the coefficient source is:

```text
StructuredTreeBlocks
```

### Serial Restrictions

The serial solver checks:

- `comm_size == 1`,
- sparse Jacobian is completed if using `SparseJacobian`,
- structured `TreeLinearization` is attached if using `StructuredTreeBlocks`,
- structured row and dof counts match metadata,
- all residual rows are locally available,
- all correction dofs are locally available.

### Symbolic Plan

The serial solver builds a symbolic plan in its constructor. This plan is reused across Newton corrections.

For each element, `ElementSolvePlan` records:

- element index and global id,
- inlet pressure local dof,
- inlet flow unknown index,
- whether the element is a leaf,
- unknown global dof ids,
- unknown local dof ids,
- equation rows used in the local dense block,
- child-interface plans,
- a text context used in error messages.

The element-local unknown vector contains every element dof except inlet pressure. Inlet pressure is treated as an external scalar parameter.

For a leaf element, the local rows are:

```text
element state-equation rows
leaf outlet boundary row
```

For an internal element, the local rows are:

```text
element state-equation rows
downstream junction flow-conservation row
```

Pressure-continuity junction rows are not part of the dense local block. They are used to express child inlet pressure as an affine function of parent outlet pressure.

The solver also preallocates per-element workspaces:

- dense local matrix,
- constant RHS,
- inlet-pressure RHS,
- intercept solution vector,
- slope solution vector,
- child pressure slopes,
- child pressure intercepts.

It stores one `SubtreeRelation` per element and one inlet-pressure value per element for top-down recovery.

### Bottom-Up Subtree Condensation

The serial solver processes `tree_metadata.bottom_up_layers` from leaves to root.

For each element, it builds a dense local system:

```text
A * u = b_const + b_pin * delta_p_in
```

where:

- `u` contains all element dof corrections except inlet pressure,
- `delta_p_in` is the element inlet pressure correction,
- `A` contains coefficients for the local unknowns,
- `b_const` contains `-F` and child-substitution contributions,
- `b_pin` is the negative coefficient of inlet pressure in each local row.

The code solves the dense system for two right-hand sides at once:

```text
A * intercept = b_const
A * slope     = b_pin
```

Thus:

```text
u = intercept + slope * delta_p_in
```

The element's subtree relation is read from the inlet-flow unknown entry:

```text
G = slope[inlet_flow_unknown_index]
h = intercept[inlet_flow_unknown_index]
```

so:

```text
delta_q_in = G * delta_p_in + h
```

### Child Subtree Substitution

For an internal element, each child subtree has already produced:

```text
delta_q_child_in = G_child * delta_p_child_in + h_child
```

The pressure-continuity junction row has the current form:

```text
a_parent * delta_p_parent_out + a_child * delta_p_child_in = rhs_pressure
```

where:

```text
rhs_pressure = -F_pressure_row
```

The solver computes:

```text
delta_p_child_in = slope_child_pressure * delta_p_parent_out + intercept_child_pressure
```

with:

```text
slope_child_pressure = -a_parent / a_child
intercept_child_pressure = rhs_pressure / a_child
```

Then it substitutes the child subtree relation into the parent flow-conservation row:

```text
delta_q_child_in
  = G_child * (slope_child_pressure * delta_p_parent_out + intercept_child_pressure) + h_child
```

This modifies the parent local dense row for `delta_p_parent_out` and shifts the RHS by the child-flow intercept contribution. This is the local Schur complement step that removes the child subtree from the parent solve.

### Local Dense Solver

The local dense solver is an internal Gaussian-elimination helper with partial pivoting. It solves the intercept and slope systems together by overwriting the local matrix and both RHS vectors.

If the best pivot in a column is below `pivot_tolerance`, the solver throws a singular or underconstrained local block error. The default pivot tolerance is:

```text
1.0e-12
```

### Root Solution

After bottom-up condensation reaches the root, the root inlet boundary condition determines the root inlet pressure correction.

The root boundary row is:

```text
boundary_coeff * delta_p_root = -F_boundary
```

Therefore:

```text
delta_p_root = -F_boundary / boundary_coeff
```

In code, `rhs_value(residual, root_boundary_row)` returns `-F_boundary`, so the implementation computes:

```text
root_inlet_pressure = rhs_value / root_boundary_coeff
```

### Top-Down Recovery

The serial solver processes `tree_metadata.top_down_layers` from root to leaves.

For each element:

```text
known delta_p_in
u = slope * delta_p_in + intercept
write inlet pressure and all local unknown corrections into delta
```

For internal elements, the recovered parent outlet pressure is used with the stored pressure-continuity relation to compute each child inlet pressure:

```text
delta_p_child_in = child_pressure_slope * delta_p_parent_out + child_pressure_intercept
```

Those child inlet pressures drive the next top-down layer.

### Construction Of The Newton Correction Vector

The serial solver writes directly into the correction vector using global dof IDs:

```text
delta.replace_global_value(global_dof_id, value)
```

After `TreeNewtonLinearSolver::solve()` returns, `NewtonSolver` applies:

```text
x_solution_ = x_solution_ + delta
```

## Distributed Tree Solver: C++ Implementation And Data Flow

The distributed solver class is:

```text
DistributedTreeNewtonLinearSolver
```

It is constructed with:

```text
DistributedTreeNewtonLinearSolverContext{
  tree_metadata,
  locally_relevant_dof_map,
  pivot_tolerance,
  profile
}
```

Runtime `NewtonTree` selects this solver when:

```text
comm_size > 1
```

The distributed solver reports:

```text
linearization_type() == StructuredTreeBlocks
```

so the custom Newton driver assembles `TreeLinearization`, not a sparse Jacobian.

### Difference From The Serial Solver

The serial solver uses local row and local dof IDs and requires the full correction system on one rank.

The distributed solver uses global row and global dof IDs in its symbolic plan. Each rank processes only elements for which:

```text
TreeElementMetadata::owner_rank == my_rank
```

It uses a rank-local coefficient provider that maps:

```text
global row id -> local residual row id
global dof id -> local locally-relevant dof id
```

Then it reads coefficients from the rank-local `TreeLinearization`.

If a rank owns an element but does not locally have a required row or locally relevant dof, the solver fails with a clear error. This preserves the current assembly ownership assumptions instead of silently gathering missing data.

### Distributed Bottom-Up Phase

The distributed solver processes bottom-up layers collectively. In each layer:

1. Each rank processes only its owned elements in that layer.
2. Same-rank child subtree relations are read from local storage.
3. Cross-rank child subtree relations are received from child owners.
4. Each owned element is condensed to:

```text
delta_q_in = G * delta_p_in + h
```

5. If the parent element is owned by another rank, the child owner sends:

```text
element_index, G, h
```

This is the true partition-boundary version of the bottom-up tree relation exchange. The current code uses `MPI_Alltoallv` per layer to exchange only the boundary relations needed by other ranks.

### Distributed Root Solution

Only the root element owner computes the root inlet pressure correction from the root inlet boundary row:

```text
delta_p_root = -F_boundary / boundary_coeff
```

Other ranks learn their needed inlet pressures during the top-down phase.

### Distributed Top-Down Phase

The distributed solver processes top-down layers collectively. In each layer:

1. A rank processes only owned elements whose inlet pressure is known.
2. The rank recovers local element corrections using the stored slope/intercept data.
3. Same-rank child inlet pressures are stored locally.
4. Cross-rank child inlet pressures are sent to child owners as:

```text
child_element_index, child_inlet_pressure
```

These messages are also exchanged with `MPI_Alltoallv` per tree layer.

### Final Correction Scatter

Element owners compute correction values, but the `delta` vector is distributed according to the correction map. The element owner may not own every correction-vector entry it computed.

The distributed solver queries correction-vector owners with:

```text
correction_map.remote_id_list(...)
```

For each computed dof correction:

- if the correction-vector owner is the current rank, it writes directly to `delta`,
- otherwise it queues a scatter message:

```text
global_dof_id, delta_value
```

After the top-down recovery, the queued correction messages are exchanged with `MPI_Alltoallv`, and receiving ranks insert the values into their local `delta` vectors.

### MPI Communication Strategy

The current distributed tree solver communicates three message categories:

| Phase | Message | Purpose |
| --- | --- | --- |
| Bottom-up | `(element_index, G, h)` | Send condensed child subtree relation to parent owner |
| Top-down | `(child_element_index, child_inlet_pressure)` | Send known inlet pressure correction to child owner |
| Final scatter | `(global_dof_id, delta_value)` | Move computed correction entries to correction-vector owner |

The distributed solver no longer performs the earlier replicated global solve strategy. It does not gather all residual entries or all structured coefficients to all ranks. It performs rank-local condensation/recovery and communicates only partition-boundary relations, partition-boundary inlet pressures, and final correction scatter values.

Current communication implementation details:

- `MPI_Alltoall` exchanges per-rank message counts,
- `MPI_Alltoallv` exchanges IDs and scalar payloads,
- relation messages use one integer and two doubles,
- scalar messages use one integer and one double,
- communication is collective per tree layer, not asynchronous point-to-point.

The profiling counters distinguish:

- `boundary_relation_message_count`,
- `boundary_pressure_message_count`,
- `correction_scatter_message_count`,
- `communication_time`,
- `communication_bytes`.

## Main Classes And Files

### Runtime Driver And Input

| File | Main role |
| --- | --- |
| `src/reduced_lung/src/4C_reduced_lung_main.cpp` | Simulation setup, time loop, runtime solver selection |
| `src/reduced_lung/src/4C_reduced_lung_input.hpp` | `ReducedLungParameters`, including `NonlinearSolverType` |
| `src/reduced_lung/src/4C_reduced_lung_input.cpp` | Valid input registration for `nonlinear_solver` |

### Shared Assembly And NOX Path

| File | Main role |
| --- | --- |
| `src/reduced_lung/src/4C_reduced_lung_helpers.hpp` | Assembly pipeline, NOX context/solver declarations, helper declarations |
| `src/reduced_lung/src/4C_reduced_lung_helpers.cpp` | Pipeline callback registration, NOX wrapper, discretization/model/map helpers |

### Custom Newton And Linear Solver Interface

| File | Main role |
| --- | --- |
| `src/reduced_lung/src/4C_reduced_lung_newton_solver.hpp/cpp` | Full-step Newton nonlinear loop |
| `src/reduced_lung/src/4C_reduced_lung_linear_solver.hpp/cpp` | `NewtonLinearSolver` interface and sparse backend |

### Tree Solver Infrastructure

| File | Main role |
| --- | --- |
| `src/reduced_lung/src/4C_reduced_lung_tree_metadata.hpp/cpp` | Directed tree metadata and validation |
| `src/reduced_lung/src/4C_reduced_lung_tree_linearization.hpp/cpp` | Structured local derivative storage |
| `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.hpp` | Serial and distributed tree solver declarations |
| `src/reduced_lung/src/4C_reduced_lung_tree_linear_solver.cpp` | Serial bottom-up/top-down tree solver |
| `src/reduced_lung/src/4C_reduced_lung_distributed_tree_linear_solver.cpp` | MPI bottom-up/top-down tree solver |

### Physics Modules Modified For Structured Blocks

| File group | Main role |
| --- | --- |
| `src/reduced_lung/src/airways/*` | Airway residuals, sparse Jacobians, tree linearization blocks, state/history updates |
| `src/reduced_lung/src/terminal_units/*` | Terminal-unit residuals, sparse Jacobians, tree linearization blocks, state/history updates |
| `src/reduced_lung/src/4C_reduced_lung_junctions.hpp/cpp` | Junction residuals, sparse Jacobians, tree linearization blocks |
| `src/reduced_lung/src/4C_reduced_lung_boundary_conditions.hpp/cpp` | Boundary residuals, sparse Jacobians, tree linearization blocks |

### Profiling And Benchmarks

| File | Main role |
| --- | --- |
| `src/reduced_lung/src/4C_reduced_lung_solver_profile.hpp` | Optional timing and counter structs |
| `src/reduced_lung/benchmark_tests/4C_reduced_lung_solver_benchmark.cpp` | Google Benchmark harness |
| `src/reduced_lung/benchmark_tests/CMakeLists.txt` | Benchmark target registration |

### Tests Added Or Extended

| File | Main role |
| --- | --- |
| `src/reduced_lung/tests/4C_reduced_lung_newton_solver_test.np2.cpp` | Direct custom Newton analytical terminal-unit test |
| `src/reduced_lung/tests/4C_reduced_lung_newton_vs_nox_test.np2.cpp` | NOX versus custom Newton comparisons |
| `src/reduced_lung/tests/4C_reduced_lung_tree_metadata_test.cpp` | Tree metadata construction and validation tests |
| `src/reduced_lung/tests/4C_reduced_lung_tree_linear_solver_test.cpp` | Serial sparse/tree correction and full workflow comparisons |
| `src/reduced_lung/tests/4C_reduced_lung_distributed_tree_linear_solver_test.np2.cpp` | Two-rank distributed tree correction comparisons |
| `src/reduced_lung/tests/4C_reduced_lung_input_pipeline_test.cpp` | Input/default solver policy coverage |

## Validation Strategy

Validation is layered so that failures can be localized to either the nonlinear loop, sparse backend, tree metadata, tree linear solver, distributed communication, or runtime input wiring.

### Existing Component Unit Tests

Existing tests continue to cover:

- input pipeline and model creation,
- helper functions and dof layout,
- airway Jacobians versus finite differences,
- terminal-unit Jacobians versus finite differences,
- junction residual/Jacobian assembly,
- boundary-condition residual/Jacobian assembly,
- original NOX single-terminal-unit behavior.

These tests guard the preserved physics and assembly behavior.

### Custom Newton Tests

`4C_reduced_lung_newton_solver_test.np2.cpp` constructs the custom Newton solver directly for a single terminal unit with Ogden elasticity and checks the analytical volume behavior over multiple time steps.

`4C_reduced_lung_newton_vs_nox_test.np2.cpp` builds independent fixtures and compares NOX and custom Newton with sparse backend for:

- single terminal unit with Ogden elasticity,
- single terminal unit with linear elasticity,
- serial rigid airways,
- rigid-airway bifurcation.

The comparisons check solution vectors, owned dofs, residual convergence, selected flow balances, and terminal-unit volume history where relevant.

### Tree Metadata Tests

`4C_reduced_lung_tree_metadata_test.cpp` checks:

- connection metadata,
- bifurcation metadata with terminal-unit leaves,
- traversal layer construction,
- element classification,
- junction dof metadata,
- boundary-side classification,
- directed cycle rejection,
- unsupported branch degree rejection,
- missing root boundary rejection.

### Serial Sparse-Versus-Tree Correction Tests

`4C_reduced_lung_tree_linear_solver_test.cpp` compares one Newton correction from:

```text
SparseNewtonLinearSolver
sparse-source TreeNewtonLinearSolver
structured-source TreeNewtonLinearSolver
```

The correction vector is compared entry-by-entry.

Covered cases include:

- single terminal unit,
- serial rigid airways,
- rigid-airway bifurcation,
- Kelvin-Voigt airways,
- nonlinear airway resistance,
- Four-element Maxwell terminal unit,
- mixed airway/terminal-unit tree,
- reuse of one structured tree solver instance across several state/time points.

### Full Workflow Comparisons

The same serial tree test file compares complete time-step workflows for:

```text
Nox
NewtonSparse
NewtonTree
```

Covered workflow cases include:

- single terminal unit,
- serial rigid airways,
- bifurcation rigid airways,
- mixed airways and terminal units.

The tests compare final solution vectors, owned dofs, residual norms, flow balances, and terminal-unit volume updates where appropriate.

### Distributed NP2/MPI Tests

`4C_reduced_lung_distributed_tree_linear_solver_test.np2.cpp` tests the distributed tree correction solver on two ranks.

It compares:

```text
SparseNewtonLinearSolver correction
DistributedTreeNewtonLinearSolver correction
```

Covered distributed cases are:

- serial airway chain,
- bifurcation airways.

Both tests assert that the tree metadata contains at least one cross-rank edge, so the boundary communication path is actually exercised.

### YAML Runtime Tests

Runtime YAML variants were added and registered in `tests/list_of_tests.cmake` for explicit `NewtonSparse` and `NewtonTree` selection. Examples include:

- `reduced_lung_terminal_unit_newton_sparse.4C.yaml`,
- `reduced_lung_terminal_unit_newton_tree.4C.yaml`,
- `reduced_lung_serial_airways_flow_newton_sparse.4C.yaml`,
- `reduced_lung_serial_airways_flow_newton_tree.4C.yaml`,
- `reduced_lung_aw_bifurcation_flow_newton_sparse.4C.yaml`,
- `reduced_lung_aw_bifurcation_flow_newton_tree.4C.yaml`,
- `reduced_lung_3_aw_2_tu_newton_sparse.4C.yaml`,
- `reduced_lung_3_aw_2_tu_newton_tree.4C.yaml`,
- `reduced_lung_3_aw_2_tu_4elemax_and_kv_newton_sparse.4C.yaml`,
- `reduced_lung_3_aw_2_tu_4elemax_and_kv_newton_tree.4C.yaml`.

The original YAML files still omit `nonlinear_solver` and therefore still use `Nox` by default.

### Verification Commands Used During Implementation

The implementation notes record successful runs of commands such as:

```text
cmake --build build/debug --target unittests_reduced_lung --parallel 4
cmake --build build/debug --target unittests_reduced_lung.np2 --parallel 4
ctest -R "^unittests_reduced_lung(\.np2)?$" --output-on-failure
ctest -R "reduced_lung_.*newton_.*\.4C\.yaml-p1$" --output-on-failure
git diff --check
```

This document creation did not rerun the test suite. The summary above is based on the current code and the recorded phase verification notes.

## Profiling And Benchmark Infrastructure

Optional profiling structs are available for:

- `NoxSolver`,
- `NewtonSolver`,
- `SparseNewtonLinearSolver`,
- tree Newton linear solvers.

Recorded quantities include:

- full NOX solve time,
- NOX residual and Jacobian evaluation counts,
- residual assembly time,
- sparse Jacobian assembly time,
- sparse matrix completion time,
- full custom Newton solve time,
- state synchronization time,
- structured tree-linearization assembly time,
- linear correction solve time,
- sparse linear solver time,
- tree total solve time,
- tree bottom-up time,
- tree top-down time,
- tree dense local solve time and count,
- tree coefficient lookup time and count,
- tree symbolic workspace counters,
- last residual norm history,
- last increment norm history,
- distributed communication time, bytes, and message counts.

The benchmark harness is registered through:

```text
four_c_auto_define_benchmark_tests()
```

Benchmark groups include:

- `ReducedLung/FullSolve/...` comparing `Nox`, `NewtonSparse`, and `NewtonTree`,
- `ReducedLung/AssemblyPhases/BalancedAirways`,
- `ReducedLung/LinearSolve/BalancedAirways/Sparse`,
- `ReducedLung/LinearSolve/BalancedAirways/StructuredTree`.

Benchmark cases include:

- single terminal unit,
- 9-element serial airway chain,
- 4-level balanced airway tree,
- balanced airway trees with 2, 3, 4, and 5 levels for linear-solve scaling.

The benchmark code currently enforces serial execution with `ensure_serial_benchmark()`. Distributed tree benchmark coverage has not yet been added.

The implementation notes record a debug-build benchmark run that showed much lower serial linear-solve times for the structured tree path than the sparse path. Those numbers are useful as a sanity check only. They are from a debug build and should not be treated as final performance data.

## Current Limitations And Bottlenecks

### Solver Policy Limitations

`Nox` remains the default. This is intentional. `NewtonSparse` and `NewtonTree` are opt-in workflows until validation and performance evidence are broader.

NOX remains a required reduced-lung module dependency because the NOX path is still selectable and still the default.

### Nonlinear Algorithm Limitations

`NewtonSolver` is full-step Newton only. It has no damping, no line search, and no trust-region globalization. This is close to the configured NOX full-step behavior, but it does not reproduce all of NOX's nonlinear solver infrastructure.

The convergence checks are simple absolute two-norm checks on residual and increment. There is no scaling or adaptive tolerance logic.

### Tree Topology Limitations

The tree metadata and tree solver require:

- one connected directed tree,
- exactly one root,
- acyclic parent-child topology,
- branch degree at most two,
- terminal units as leaves,
- root inlet boundary closure,
- outlet boundary closure for every leaf,
- square equation/dof count.

Unsupported layouts throw. No automatic fallback to `NewtonSparse` is implemented for explicit `NewtonTree` runs.

Element direction is critical. The solver assumes:

```text
element_nodes[0] = inlet / parent side
element_nodes[1] = outlet / child side
```

### Data-Layout Limitations

The current duplicated endpoint dof formulation is preserved. No shared-node pressure/flow formulation has been introduced.

The nonlinear solution vector `x_` is still created on the row map, while owned dofs use the locally owned dof map. The implementation relies on the square system and compatible global numbering so exports and correction updates work consistently.

### Serial Tree Performance Bottlenecks

The serial tree solver already reuses a symbolic solve plan and per-element dense workspaces. However, runtime `NewtonTree` still rebuilds the numeric `TreeLinearization` every Newton correction.

Coefficient access is still through row/dof lookup in `TreeLinearization`. This is much lighter than sparse matrix assembly and sparse factorization, but it is not yet a fully preplanned dense block update path.

The local dense solver is simple Gaussian elimination. Current local blocks are small, but this has not been replaced with specialized closed-form kernels or batched dense kernels.

### Distributed Tree Bottlenecks

The distributed solver now communicates only partition-boundary data, but it still uses collective `MPI_Alltoallv` exchanges per tree layer.

Current bottlenecks and limitations include:

- no asynchronous point-to-point communication,
- no precomputed per-layer communication plan,
- correction owner lookup is queried during each solve,
- assumptions that owned element rows and required locally relevant dofs are locally available,
- distributed validation currently covers correction solves for two airway topologies on two ranks,
- distributed mixed airway/terminal-unit workflows are not yet broadly validated,
- no distributed benchmark harness exists yet.

### Benchmark Limitations

The recorded benchmark numbers are from a debug build. Release-build measurements are still needed for reliable performance conclusions.

External heap allocation profiling was recommended but not recorded in the implementation notes.

## Logical Next Optimization Steps

The following items are proposed future work. They are not already implemented behavior.

### Broaden Validation Before Changing Defaults

Recommended validation work:

- run `NewtonTree` on more representative full YAML inputs,
- add distributed full-workflow tests, not only distributed correction tests,
- add distributed mixed airway/terminal-unit cases,
- add larger realistic tree cases,
- compare `Nox`, `NewtonSparse`, and `NewtonTree` over multiple time steps and model combinations,
- track nonlinear iteration counts, residual histories, and increment histories,
- keep `Nox` as default until these comparisons are stable.

### Improve Structured Block Assembly

Current runtime `NewtonTree` avoids sparse Jacobian assembly, but still fills a generic row-oriented `TreeLinearization`.

A faster path would precompute a symbolic block plan and assemble directly into per-element, per-junction, and per-boundary dense blocks used by the tree solver.

Possible steps:

- keep derivative formulas in the physics modules,
- replace generic row/dof `set_value` calls with structured block writers,
- store direct indices into tree-solver workspaces,
- reuse numeric block storage across Newton iterations,
- avoid repeated coefficient lookup by row and dof.

### Optimize Distributed Communication

The distributed algorithm is now partition-boundary based, but the communication implementation can be improved.

Possible steps:

- cache correction owner lookup in the distributed solver,
- precompute per-layer send and receive plans from tree metadata,
- replace all-to-all layer exchanges with point-to-point communication between actual neighboring tree partitions,
- overlap local condensation/recovery with communication where possible,
- add NP4 and larger distributed benchmarks.

### Improve Nonlinear Robustness If Needed

If full-step Newton is not robust enough on larger or more nonlinear cases, add controlled globalization to the custom Newton path.

Possible options:

- damping factor,
- backtracking line search,
- residual-based step acceptance,
- fallback from `NewtonTree` to `NewtonSparse` only if explicitly requested by solver policy.

### Produce Final Performance Data

Recommended measurement work:

- run release-build benchmarks,
- separate residual assembly, structured block assembly, tree solve, and communication costs,
- measure sparse factorization separately from sparse assembly,
- profile heap allocations with external tools,
- benchmark serial and distributed scaling on increasing tree sizes,
- report speedups together with nonlinear iteration counts and residual histories.

### Default Solver Policy

Do not make `NewtonTree` the default yet. The current policy should remain:

```text
Nox          = default/reference
NewtonSparse = explicit sparse custom-Newton reference/debug path
NewtonTree   = explicit structured tree-solver path
```

Changing the default should require broad validation, release-build performance evidence, and a clear fallback or support policy for unsupported tree layouts and MPI distributions.
