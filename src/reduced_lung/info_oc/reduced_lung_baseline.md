# reduced_lung Baseline Documentation

This document is a baseline snapshot of the current `src/reduced_lung/` implementation before the thesis changes. It describes the current NOX-based/generic solver workflow as it exists now. It is not intended to describe the final future implementation. Future work will likely replace or modify parts of the current NOX-based solver workflow with a custom tree-based solver, but this document does not specify or implement that future solver.

This document intentionally ignores `src/1d_pipe_flow/`, which is a separate reduced-lung 1D pipe-flow path.

## Purpose

The `reduced_lung` module implements a reduced-dimensional lung tree solver. The current implementation builds a lightweight 1D tree discretization from input topology, instantiates per-element airway and terminal-unit models, adds node-level coupling and boundary equations, and solves the resulting nonlinear algebraic system at each time step using the 4C NOX nonlinear solver infrastructure.

The modeled tree consists of line elements with two endpoint nodes. Each element is either an airway or a terminal unit. Airways carry pressure and flow unknowns and may use rigid or Kelvin-Voigt wall mechanics. Terminal units carry pressure and flow unknowns and use a rheological model plus an elasticity model. Internal tree nodes create connection or bifurcation equations. Boundary nodes create pressure or flow boundary-condition equations.

## Main Entry Points

The main reduced-lung entry point is `ReducedLung::reduced_lung_main()` in `src/4C_reduced_lung_main.cpp` and `src/4C_reduced_lung_main.hpp`.

Important external wiring:

- `apps/global_full/4C_global_full_entrypoint_switch.cpp` calls `ReducedLung::reduced_lung_main()` for `Core::ProblemType::reduced_lung`.
- `src/global_legacy_module/4C_global_legacy_module_problem_type_string.hpp` maps input problem type string `Reduced_Lung` to `Core::ProblemType::reduced_lung`.
- `src/global_legacy_module/4C_global_legacy_module_validparameters.cpp` registers `ReducedLung::valid_parameters()` with the global valid-parameter list.
- `src/global_data/4C_global_data_read.cpp` creates an empty `red_airway` discretization for `Core::ProblemType::reduced_lung`, but the current reduced-lung solver also builds its own local `Core::FE::Discretization` named `reduced_lung` from the `reduced_dimensional_lung` topology.

## Directory Layout

- `CMakeLists.txt`: defines the reduced-lung module and declares internal dependencies on `core`, `global_data`, and `solver_nonlin_nox`.
- `src/CMakeLists.txt`, `src/airways/CMakeLists.txt`, `src/terminal_units/CMakeLists.txt`: use `four_c_auto_define_module()`.
- `src/4C_reduced_lung_main.cpp`, `src/4C_reduced_lung_main.hpp`: top-level simulation orchestration.
- `src/4C_reduced_lung_input.cpp`, `src/4C_reduced_lung_input.hpp`: input data structures and valid-parameter specification for `reduced_dimensional_lung`.
- `src/4C_reduced_lung_helpers.cpp`, `src/4C_reduced_lung_helpers.hpp`: discretization construction, model creation helpers, global/local maps, NOX wrapper, output collection, and assembly pipeline.
- `src/4C_reduced_lung_boundary_conditions.cpp`, `src/4C_reduced_lung_boundary_conditions.hpp`: boundary-condition model creation and assembly.
- `src/4C_reduced_lung_junctions.cpp`, `src/4C_reduced_lung_junctions.hpp`: connection and bifurcation equations at internal nodes.
- `src/airways/`: airway model data, registries, flow resistance, wall mechanics, residual/Jacobian callbacks, and end-of-step updates.
- `src/terminal_units/`: terminal-unit model data, registries, rheology, elasticity, residual/Jacobian callbacks, and end-of-step updates.
- `src/4C_reduced_lung_aaa_approx.cpp`, `src/4C_reduced_lung_aaa_approx.hpp`: standalone AAA rational approximation utility. It is tested in this module but is not part of the main NOX time-stepping workflow in `ReducedLungSimulation`.
- `tests/`: unit and integration-style tests for input, helpers, airways, terminal units, boundary conditions, junctions, AAA, and the NOX solver.

## Input Model

Input is stored in `ReducedLung::ReducedLungParameters` from `src/4C_reduced_lung_input.hpp`. `ReducedLung::valid_parameters()` in `src/4C_reduced_lung_input.cpp` registers the `reduced_dimensional_lung` section.

The main input groups are:

- `dynamics`: `time_increment`, `number_of_steps`, `restart_every`, `results_every`, `linear_solver`, `max_nonlinear_iterations`, `nonlinear_residual_tolerance`, `nonlinear_increment_tolerance`, and `output_verbosity`.
- `lung_tree.topology`: `num_nodes`, `num_elements`, `node_coordinates`, and `element_nodes`.
- `lung_tree.element_type`: per-element `Airway` or `TerminalUnit`.
- `lung_tree.generation`: per-element generation value. It is parsed and tested but is not central in the current solver assembly.
- `lung_tree.airways`: airway radius, flow resistance model, inertia flag, wall model, and wall parameters.
- `lung_tree.terminal_units`: terminal-unit rheology and elasticity model selections and parameters.
- `boundary_conditions`: number of conditions, type, node id, and either function id or constant value.
- `air_properties`: air density and dynamic viscosity.

Representative full input files exist outside this directory in `tests/input_files/`, for example `reduced_lung_3_aw_2_tu.4C.yaml`, `reduced_lung_aw_bifurcation_flow.4C.yaml`, `reduced_lung_serial_airways_flow.4C.yaml`, and `reduced_lung_terminal_unit.4C.yaml`. These use `PROBLEMTYPE: "Reduced_Lung"`, a `reduced_dimensional_lung` section, and a `SOLVER n` section referenced by `dynamics.linear_solver`.

## Initialization Workflow

The current simulation is implemented by the private `ReducedLungSimulation` class in `src/4C_reduced_lung_main.cpp`.

`ReducedLung::reduced_lung_main(Global::Problem& problem)` performs these high-level steps:

- Calls `make_reduced_lung_context_from_problem(problem)` to gather `ReducedLungParameters`, MPI communicator, rebalance parameters, IO parameters, linear solver parameters, solver parameter callback, output control, and function manager.
- Constructs `ReducedLungSimulation`.
- Calls `ReducedLungSimulation::initialize()`.
- Calls `ReducedLungSimulation::run()`.

`ReducedLungSimulation::initialize()` performs this sequence:

- `validate_parameters()` checks positive `time_increment`, non-negative `number_of_steps`, positive `results_every`, and positive `max_nonlinear_iterations`.
- `build_discretization()` calls `build_discretization_from_topology()`, then `actdis_->fill_complete()`, then creates `Core::IO::DiscretizationVisualizationWriterMesh` for runtime VTK output.
- `build_element_models()` creates airway and terminal-unit model containers, creates global dof offsets, assigns global dof ids, and creates evaluator callbacks.
- `build_node_entities()` builds node-to-element adjacency, creates boundary-condition models/evaluators, creates junction data, and prints instantiated object counts.
- `assign_equation_ids()` assigns local row/equation ids to airways, terminal units, junctions, and boundary conditions.
- `build_maps_and_local_ids()` creates the domain, row, and column maps, assigns global equation ids to junctions and boundary conditions, and assigns local dof ids to all model blocks.
- `build_linear_system_and_solver()` creates vectors/matrix, creates the default NOX assembly pipeline, builds `NoxSolverContext`, and constructs `NoxSolver`.

## Discretization and Topology

`build_discretization_from_topology()` in `src/4C_reduced_lung_helpers.cpp` builds a minimal 3D `Core::FE::Discretization` from input topology:

- Only MPI rank 0 adds nodes and elements to `Core::FE::DiscretizationBuilder<3>`.
- Input node ids and element ids are treated as 1-based in the input fields and converted to 0-based internal/global ids.
- Each element is a `Core::FE::CellType::line2` element with two nodes.
- The function validates positive node/element counts, 3-component node coordinates, 2-node element connectivity, valid node id ranges, and no self-connected elements.
- The builder distributes/rebalances the discretization using `Core::Rebalance::RebalanceParameters`.

Important convention: element node order is semantically meaningful. `element_nodes[0]` is treated as inlet/parent side and `element_nodes[1]` as outlet/child side. Junction creation and boundary dof selection depend on this order.

## Element Model Creation

`create_local_element_models()` in `src/4C_reduced_lung_helpers.cpp` loops over `discretization.my_row_element_range()` and instantiates locally owned model data.

For airway elements:

- Reads `parameters.lung_tree.element_type` and dispatches when it is `ReducedLungParameters::LungTree::ElementType::Airway`.
- Reads airway `resistance_type` and `wall_model_type`.
- Calls `Airways::ModelRegistry::add_airway_with_model_selection()`.
- Stores the per-element dof count as `2 + n_state_equations`.
- Rigid-wall airway models have one state equation and three dofs: `p1`, `p2`, `q1`.
- Kelvin-Voigt airway models have two state equations and four dofs: `p1`, `p2`, `q1`, `q2`.

For terminal-unit elements:

- Reads terminal-unit `rheological_model_type` and `elasticity_model_type`.
- Calls `TerminalUnits::ModelRegistry::add_terminal_unit_with_model_selection()`.
- Stores the per-element dof count as `3`.
- Terminal-unit dofs are `p1`, `p2`, and `q`.

After local element model creation:

- `create_global_dof_maps()` reduces local dof counts to a global `global_dof_per_ele` map and creates `first_global_dof_of_ele` by accumulating in sorted element-id order.
- `assign_global_dof_ids_to_models()` fills airway `gid_p1`, `gid_p2`, `gid_q1`, optional `gid_q2`, and terminal-unit `gid_p1`, `gid_p2`, `gid_q`.
- `TerminalUnits::create_evaluators()` and `Airways::create_evaluators()` attach concrete residual, Jacobian, internal-state, end-of-timestep, and output callbacks.

## Airway Model Code

The airway code is under `src/airways/`.

Important types:

- `Airways::AirwayData` in `4C_reduced_lung_airways_common.hpp`: struct-of-arrays storage for element ids, row ids, global/local dof ids, reference length/area, air properties, previous-step pressures/flows, and `n_state_equations`.
- `Airways::AirwayModel` in `4C_reduced_lung_airways.hpp`: one homogeneous model block containing `AirwayData`, `FlowModel`, `WallModel`, and evaluator callbacks.
- `Airways::AirwayContainer`: local vector of airway model blocks.
- `Airways::FlowModel`: `std::variant<LinearResistive, NonLinearResistive>`.
- `Airways::WallModel`: `std::variant<RigidWall, KelvinVoigtWall>`.

Model registry:

- `Airways::ModelRegistry::add_airway_with_model_selection()` in `4C_reduced_lung_airways_model_registry.cpp` selects a factory from a registry keyed by `(FlowModelType, WallModelType)`.
- Supported pairs are `Linear/Rigid`, `Linear/KelvinVoigt`, `NonLinear/Rigid`, and `NonLinear/KelvinVoigt`.
- `register_or_access_airway_model()` groups elements sharing the same flow/wall model pair into one `AirwayModel` block.
- `add_airway_element()` computes element length from node coordinates, reference area from radius, stores air properties, initializes previous-step states to zero, and appends flow/wall parameters.

Flow resistance code:

- `4C_reduced_lung_airways_flow_resistance.hpp/cpp` defines `LinearResistive`, `NonLinearResistive`, and `FlowResistance` callback factories.
- Linear resistance uses Poiseuille resistance from `ComputePoiseuilleResistance`.
- Nonlinear resistance stores `turbulence_factor_gamma` and `k_turb`, with `k_turb` updated from the current flow in `make_internal_state_updater()`.
- Inertia is optional per element through `has_inertia` and contributes `density * ref_length / area` when enabled.

Wall mechanics code:

- `4C_reduced_lung_airways_wall_mechanics.hpp/cpp` defines `RigidWall` and `KelvinVoigtWall`.
- `WallMechanics::make_residual_evaluator()` composes wall mechanics with the selected flow resistance and inertia evaluators.
- `WallMechanics::make_jacobian_evaluator()` composes analytic derivative callbacks for the selected flow/wall pair.
- Rigid walls assemble one momentum-style equation per element.
- Kelvin-Voigt walls assemble a momentum equation and a mass/wall equation per element, update wall area and derived wall quantities, and carry `area_n` history.

Top-level airway functions in `4C_reduced_lung_airways.cpp`:

- `Airways::update_residual_vector()` loops over model blocks and calls `model.residual_evaluator()`.
- `Airways::update_jacobian()` loops over model blocks and calls `model.jacobian_evaluator()`.
- `Airways::assign_local_equation_ids()` assigns local row ids and advances by `n_state_equations` per element.
- `Airways::assign_local_dof_ids()` maps global dofs to local ids in the locally relevant dof map.
- `Airways::update_internal_state_vectors()` synchronizes internal model states during nonlinear iterations.
- `Airways::end_of_timestep_routine()` advances previous-step `p1_n`, `p2_n`, `q1_n`, optional `q2_n`, and model-specific history.

## Terminal-Unit Model Code

The terminal-unit code is under `src/terminal_units/`.

Important types:

- `TerminalUnits::TerminalUnitData` in `4C_reduced_lung_terminal_unit_common.hpp`: struct-of-arrays storage for element ids, row ids, global/local dof ids, current volume `volume_v`, and reference volume `reference_volume_v0`.
- `TerminalUnits::TerminalUnitModel` in `4C_reduced_lung_terminal_unit.hpp`: one homogeneous model block containing `TerminalUnitData`, `RheologicalModel`, `ElasticityModel`, and evaluator callbacks.
- `TerminalUnits::TerminalUnitContainer`: local vector of terminal-unit model blocks.
- `TerminalUnits::RheologicalModel`: `std::variant<KelvinVoigt, FourElementMaxwell>`.
- `TerminalUnits::ElasticityModel`: `std::variant<LinearElasticity, OgdenHyperelasticity>`.

Model registry:

- `TerminalUnits::ModelRegistry::add_terminal_unit_with_model_selection()` in `4C_reduced_lung_terminal_unit_model_registry.cpp` selects a factory from a registry keyed by `(RheologicalModelType, ElasticityModelType)`.
- Supported pairs are `KelvinVoigt/Linear`, `KelvinVoigt/Ogden`, `FourElementMaxwell/Linear`, and `FourElementMaxwell/Ogden`.
- `register_or_access_terminal_unit_model()` groups elements sharing the same rheology/elasticity pair into one `TerminalUnitModel` block.
- `add_terminal_unit_element()` computes a reference/current terminal-unit volume from the element endpoint distance treated as a radius: `(4/3) * pi * r^3`.

Elasticity code:

- `4C_reduced_lung_terminal_unit_elasticity.hpp/cpp` defines `LinearElasticity` and `OgdenHyperelasticity`.
- `Elasticity::make_elastic_pressure_evaluator()` returns a callback that computes elastic pressure for the current trial dofs and `dt`.
- `Elasticity::make_elastic_pressure_gradient_evaluator()` returns a callback for `d(p_el)/dq`.
- `Elasticity::make_output_evaluator()` emits `elastic_pressure` for high output verbosity.

Rheology code:

- `4C_reduced_lung_terminal_unit_rheology.hpp/cpp` defines `KelvinVoigt` and `FourElementMaxwell`.
- `Rheology::make_residual_evaluator()` combines rheology with the selected elastic pressure callback.
- `Rheology::make_jacobian_evaluator()` combines rheology with the selected elastic pressure-gradient callback.
- Kelvin-Voigt terminal units assemble one residual equation involving `p1 - p2`, elastic pressure, viscosity, and flow.
- Four-element Maxwell terminal units add Maxwell branch parameters and a `maxwell_pressure_p_m` history variable.
- `Rheology::make_end_of_timestep_routine()` updates Maxwell history for `FourElementMaxwell`; Kelvin-Voigt has no additional rheology history update.

Top-level terminal-unit functions in `4C_reduced_lung_terminal_unit.cpp`:

- `TerminalUnits::update_residual_vector()` loops over model blocks and calls `model.residual_evaluator()`.
- `TerminalUnits::update_jacobian()` loops over model blocks and calls `model.jacobian_evaluator()`.
- `TerminalUnits::assign_local_equation_ids()` assigns one local row per terminal-unit element.
- `TerminalUnits::assign_local_dof_ids()` maps `gid_p1`, `gid_p2`, and `gid_q` to local ids.
- `TerminalUnits::update_internal_state_vectors()` calls model internal-state callbacks.
- `TerminalUnits::end_of_timestep_routine()` advances `volume_v += q * dt`, then calls the rheology-specific end-of-timestep callback.
- `append_volume_output()` emits `volume` for medium or higher output verbosity.

## Junction Code

Junction logic is in `src/4C_reduced_lung_junctions.cpp` and `src/4C_reduced_lung_junctions.hpp`.

Important types:

- `Junctions::ConnectionData`: stores two-element connections, equation ids, element ids, and dof ids for `{p_out_parent, p_in_child, q_out_parent, q_in_child}`.
- `Junctions::BifurcationData`: stores one-parent/two-child bifurcations, equation ids, element ids, and dof ids for `{p_out_parent, p_in_child_1, p_in_child_2, q_out_parent, q_in_child_1, q_in_child_2}`.

`Junctions::create_junctions()` loops over local row elements and inspects the outlet node, `nodes[1]`, using the global node-to-element adjacency map.

- If outlet adjacency size is 2, it creates a connection.
- If outlet adjacency size is 3, it creates a bifurcation.
- If outlet adjacency size is greater than 3, it throws.
- It assumes the first element in the adjacency list is the parent and subsequent elements are children. The correctness of this depends on topology/node ordering.

Junction equation assembly:

- `Junctions::assign_junction_local_equation_ids()` assigns two rows per connection and three rows per bifurcation.
- `Junctions::assign_junction_global_equation_ids()` obtains global row ids from the row map.
- `Junctions::assign_junction_local_dof_ids()` maps global dofs to local column ids.
- `Junctions::update_residual_vector()` assembles pressure continuity and flow conservation equations.
- `Junctions::update_jacobian()` inserts constant coupling Jacobian entries. It returns immediately if the matrix is already filled.

Connection residuals are:

- `p_out_parent - p_in_child`.
- `q_out_parent - q_in_child`.

Bifurcation residuals are:

- `p_out_parent - p_in_child_1`.
- `p_out_parent - p_in_child_2`.
- `q_out_parent - q_in_child_1 - q_in_child_2`.

## Boundary-Condition Code

Boundary-condition logic is in `src/4C_reduced_lung_boundary_conditions.cpp` and `src/4C_reduced_lung_boundary_conditions.hpp`.

Important types:

- `BoundaryConditions::BoundaryConditionData`: struct-of-arrays storage for boundary node ids, attached element ids, input ids, local/global equation ids, and constrained global/local dof ids.
- `BoundaryConditions::BoundaryConditionModel`: homogeneous group by boundary type, value source, and function id. It stores optional constant values or a cached `Core::Utils::FunctionOfTime` pointer.
- `BoundaryConditions::BoundaryConditionContainer`: vector of boundary-condition model blocks.

`BoundaryConditions::create_boundary_conditions()`:

- Validates non-negative `num_conditions`.
- Converts input `bc_node_id` from 1-based to 0-based.
- Requires the boundary node to exist in the global adjacency map and to be adjacent to exactly one element.
- Only owns a boundary condition on the rank that owns the attached element.
- Detects whether the node is the inlet or outlet endpoint of the attached element.
- Maps pressure boundary conditions to the inlet/outlet pressure dof offset, `0` or `1`.
- Maps inlet flow boundary conditions to dof offset `2`.
- Maps outlet flow boundary conditions to the last dof of the element.
- Allows at most one boundary condition per `(node, type)`.
- Groups constant boundary values together by type and groups function boundary values by type and function id.

Boundary-condition assembly:

- `BoundaryConditions::assign_local_equation_ids()` assigns one row per boundary condition.
- `BoundaryConditions::assign_global_equation_ids()` obtains global row ids from the row map.
- `BoundaryConditions::assign_local_dof_ids()` maps constrained global dof ids to local column ids.
- `BoundaryConditions::create_evaluators()` selects constant-value or function-value residual evaluators and a diagonal Jacobian evaluator.
- `BoundaryConditions::update_residual_vector()` assembles `dof_value - prescribed_value`.
- `BoundaryConditions::update_jacobian()` inserts a single diagonal value `1.0` per boundary equation and returns immediately if the matrix is already filled.

## Maps, Vectors, and Matrix Layout

Map creation is in `src/4C_reduced_lung_helpers.cpp`.

`create_domain_map()` creates the locally owned dof map. It includes:

- Airway `gid_p1`, `gid_p2`, `gid_q1`, and optional `gid_q2`.
- Terminal-unit `gid_p1`, `gid_p2`, and `gid_q`.

`create_row_map()` creates the locally owned equation-row map. It concatenates three blocks:

- State equations for local airway and terminal-unit elements.
- Coupling equations for local connections and bifurcations.
- Boundary equations for local boundary conditions.

`create_column_map()` creates the locally relevant dof/column map used by local residual and Jacobian assembly. It includes:

- Local element dofs.
- Dofs needed by local connection equations.
- Dofs needed by local bifurcation equations.
- Dofs constrained by local boundary conditions.
- It sorts and removes duplicate global dof ids.

`ReducedLungSimulation::build_linear_system_and_solver()` creates:

- `dofs_`: `Core::LinAlg::Vector<double>` on `locally_owned_dof_map_`.
- `locally_relevant_dofs_`: ghosted/relevant `Core::LinAlg::Vector<double>` on `locally_relevant_dof_map_`.
- `x_`: NOX solution vector on `row_map_`.
- `sysmat_`: `Core::LinAlg::SparseMatrix` with row map `row_map_`, column map `locally_relevant_dof_map_`, and estimated 3 entries per row.

Important caveat: the NOX solution vector `x_` is built on `row_map_`, while the owned dof vector is built on `locally_owned_dof_map_`. The current implementation assumes the equation count/order and unknown count/order are compatible enough for `Core::LinAlg::export_to(x, dofs_)`. Future agents should be careful around this if changing equation or dof ordering.

## Residual and Jacobian Assembly

The current assembly path is callback-based through `NoxAssemblyPipeline` in `src/4C_reduced_lung_helpers.hpp`.

`create_default_nox_assembly_pipeline()` registers callbacks in this order:

Residual assemblers:

- `Airways::update_residual_vector()`.
- `TerminalUnits::update_residual_vector()`.
- `Junctions::update_residual_vector()`.
- `BoundaryConditions::update_residual_vector()`.

Jacobian assemblers:

- `Airways::update_jacobian()`.
- `TerminalUnits::update_jacobian()`.
- `Junctions::update_jacobian()`.
- `BoundaryConditions::update_jacobian()`.

State updaters:

- `Airways::update_internal_state_vectors()`.
- `TerminalUnits::update_internal_state_vectors()`.

The main concrete assembly locations are:

- Airway residuals and Jacobians: `src/airways/4C_reduced_lung_airways_wall_mechanics.cpp`, with flow-resistance derivatives in `src/airways/4C_reduced_lung_airways_flow_resistance.cpp`.
- Terminal-unit residuals and Jacobians: `src/terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp`, with elastic pressure and gradients in `src/terminal_units/4C_reduced_lung_terminal_unit_elasticity.cpp`.
- Junction residuals and Jacobians: `src/4C_reduced_lung_junctions.cpp`.
- Boundary residuals and Jacobians: `src/4C_reduced_lung_boundary_conditions.cpp`.
- Runtime output vectors: `collect_runtime_output_data()` in `src/4C_reduced_lung_helpers.cpp`, plus model-specific output callbacks in airway and terminal-unit files.

Matrix fill convention:

- Some Jacobian code inserts values when `!target.filled()` and replaces nonlinear entries when the matrix is already filled.
- Junction and boundary-condition Jacobians are constant and return without work if `sysmat.filled()`.
- `NoxSolver::jacobian()` calls each Jacobian assembler and then calls `jac_matrix->complete()` if the sparse matrix is not already filled.
- Because the same `sysmat_` object is reused by NOX, future changes must respect the insert-vs-replace behavior and sparsity pattern assumptions.

## Current Solver Workflow

`ReducedLungSimulation::run()` performs time integration:

- It prints a start banner on rank 0.
- It loops from `step = 1` to `n_timesteps_`.
- For each step it calls `solve_timestep(step)`.
- It calls `write_output_if_due(step)` when `step % results_every == 0`.

`ReducedLungSimulation::solve_timestep()`:

- Prints the current timestep on rank 0.
- Advances `current_time_ += dt_`.
- Calls `nox_solver_->solve(current_time_)`.
- Calls `TerminalUnits::end_of_timestep_routine(terminal_units_, *locally_relevant_dofs_, dt_)`.
- Calls `Airways::end_of_timestep_routine(airways_, *locally_relevant_dofs_, dt_)`.

`ReducedLungSimulation::write_output_if_due()`:

- Resets `visualization_writer_`.
- Calls `collect_runtime_output_data()`.
- Writes to disk with `visualization_writer_->write_to_disk(current_time_, step)`.

The current workflow has no restart writing/reading implemented in this file, even though `restart_every` exists in the input parameters.

## NOX Usage

The reduced-lung NOX wrapper is `ReducedLung::NoxSolver` in `src/4C_reduced_lung_helpers.hpp/cpp`.

`NoxSolver` owns or references:

- `x_solution_`: NOX solution vector from `NoxSolverContext::x`.
- `dofs_`: owned dof vector.
- `locally_relevant_dofs_`: ghosted/relevant dof vector used by assemblers.
- `assembly_pipeline_`: residual/Jacobian/state callbacks.
- `dt_` and `current_time_`.
- `linear_solver_`: `std::shared_ptr<Core::LinAlg::Solver>`.
- `adapter_`: `std::optional<NOX::Nln::Adapter>`.

`create_nox_parameter_list()` configures NOX as follows:

- `Nonlinear Solver = Line Search Based`.
- Direction method `Newton`.
- Line search method `Full Step`.
- Status test is an OR of convergence and maximum iterations.
- Convergence is an AND of absolute `NormF` and absolute `NormUpdateSkipFirstIter`.
- Tolerances come from `dynamics.nonlinear_residual_tolerance` and `dynamics.nonlinear_increment_tolerance`.
- Maximum iterations comes from `dynamics.max_nonlinear_iterations`.

`NoxSolver::NoxSolver()`:

- Constructs `Core::LinAlg::Solver` from the selected linear solver parameter list and solver parameter callback.
- Builds NOX residual and Jacobian callbacks bound to `NoxSolver::residual()` and `NoxSolver::jacobian()`.
- Creates a solver map with `NOX::Nln::sol_generic` mapped to the linear solver.
- Constructs `NOX::Nln::Adapter` with the communicator, NOX parameters, solver map, `x`, Jacobian operator, and callbacks.

`NoxSolver::solve(double time)`:

- Sets `current_time_`.
- Calls `adapter_->solve()`.
- Calls `sync_state_from_x(x_solution_)` after convergence.
- Returns the NOX nonlinear iteration count.

`NoxSolver::residual()`:

- Calls `sync_state_from_x(x)` so assemblers see the current trial solution.
- Runs every residual callback in `assembly_pipeline_.residual_assemblers`.
- Returns `true` to NOX.

`NoxSolver::jacobian()`:

- Calls `sync_state_from_x(x)`.
- Requires the operator to be a `Core::LinAlg::SparseMatrix` by `dynamic_cast`.
- Runs every Jacobian callback in `assembly_pipeline_.jacobian_assemblers`.
- Completes the sparse matrix if needed.
- Returns `true` to NOX.

`NoxSolver::sync_state_from_x()`:

- Exports NOX vector `x` to `dofs_`.
- Exports `dofs_` to `locally_relevant_dofs_`.
- Runs every state updater callback in `assembly_pipeline_.state_updaters`.

External NOX adapter details:

- `NOX::Nln::Adapter` is defined in `src/solver_nonlin_nox/4C_solver_nonlin_nox_adapter.hpp/cpp`.
- It wraps a residual callback and Jacobian callback in an internal `AdapterInterface` implementing NOX required/Jacobian interfaces.
- Its `solve()` resets the NOX solver with the bound vector, calls `nox_solver_->solve()`, checks final status, copies the final NOX solution back into the bound 4C vector, and returns the iteration count.

## Linear Solver Connection

The reduced-lung linear solver is connected through the standard 4C `Core::LinAlg::Solver` infrastructure.

Parameter flow:

- Input `reduced_dimensional_lung.dynamics.linear_solver` is an integer solver id.
- `make_reduced_lung_context_from_problem()` calls `problem.solver_params(parameters.dynamics.linear_solver)` and stores the resulting `Teuchos::ParameterList` in `ReducedLungContext::linear_solver_parameters`.
- It also stores `problem.solver_params_callback()` for nested or id-based solver parameter lookup.
- `NoxSolver::NoxSolver()` constructs `Core::LinAlg::Solver(context.linear_solver_parameters, context.comm, context.solver_params_callback, Core::IO::Verbositylevel::minimal)`.
- That solver is passed to NOX in a `std::map<NOX::Nln::SolutionType, Teuchos::RCP<Core::LinAlg::Solver>>` under key `NOX::Nln::sol_generic`.

`Core::LinAlg::Solver` supports direct solvers such as `KLU2`, `MUMPS`, `UMFPACK`, `Superlu`, and iterative `Belos`, depending on the `SOLVER` setting. The reduced-lung regression inputs commonly use `SOLVER 1` with `SOLVER: "UMFPACK"` and `NAME: "Reduced_lung_solver"`.

## Outputs

Output is collected by `collect_runtime_output_data()` in `src/4C_reduced_lung_helpers.cpp` and written by `Core::IO::DiscretizationVisualizationWriterMesh`.

Always-written element fields:

- `p_1`: inlet/start pressure dof.
- `p_2`: outlet/end pressure dof.
- `q_in`: inlet/start flow dof.
- `q_out`: outlet/end flow dof. For rigid airway elements, `q_out` uses `q1`; for Kelvin-Voigt airway elements, it uses `q2`. Terminal units currently set `q_in`, while `q_out` remains `NaN` unless another model writes it.

Additional fields depend on `ReducedLungParameters::OutputVerbosity`:

- `minimal`: only `p_1`, `p_2`, `q_in`, `q_out`.
- `medium`: includes model outputs such as airway `area` and terminal-unit `volume`.
- `high`: includes diagnostics such as nonlinear airway `flow_k_turb`, terminal-unit `elastic_pressure`, and Four-element Maxwell `maxwell_pressure`.

The helper `RuntimeOutputCollector` creates element-row-map vectors initialized to `NaN` and lets model-specific output callbacks fill values for elements that support a field.

## AAA Utility

`src/4C_reduced_lung_aaa_approx.hpp/cpp` implements a standalone Adaptive Antoulas-Anderson rational approximation utility:

- `AAAOptions`: tolerance and maximum number of support points.
- `AAAResult`: support points `z`, values `f`, weights `w`, error history `errvec`, evaluation operator, pole/zero computation, and gain computation.
- `aaa()`: builds a barycentric rational approximant from real sample points and a target function callback.
- `compute_residues()`: computes residues at poles using barycentric data.

This utility has extensive tests in `tests/4C_reduced_lung_aaa_approx_test.cpp`. It is part of this module but not used by `ReducedLungSimulation` in the current baseline solver workflow.

## Tests and Build Notes

The module is auto-defined through CMake:

- `src/reduced_lung/CMakeLists.txt` calls `four_c_auto_define_module()` and declares dependencies on `core`, `global_data`, and `solver_nonlin_nox`.
- `src/reduced_lung/tests/CMakeLists.txt` calls `four_c_auto_define_tests()`.

Relevant tests in `src/reduced_lung/tests/`:

- `4C_reduced_lung_input_pipeline_test.cpp`: builds discretization and model containers from representative parameters.
- `4C_reduced_lung_helpers_test.cpp`: tests dof maps, local/global ids, local model creation, and node adjacency.
- `4C_reduced_lung_airways_test.cpp`: checks airway Jacobians against finite differences and tests airway model registry behavior.
- `4C_reduced_lung_terminal_unit_test.cpp`: checks terminal-unit Jacobians against finite differences for model combinations and tests registry errors.
- `4C_reduced_lung_boundary_conditions_test.cpp`: tests boundary condition creation, residual/Jacobian assembly, function values, and validation errors. Some tests skip when not running serially.
- `4C_reduced_lung_junctions_test.cpp`: tests connection/bifurcation residuals, Jacobians, id assignment, creation, and validation errors. Some tests skip when not running serially.
- `4C_reduced_lung_nox_solver_test.np2.cpp`: integration-style NOX solver test for a single terminal unit with Ogden elasticity. The `.np2` suffix indicates it is intended for a two-process test setup in the 4C test conventions.
- `4C_reduced_lung_aaa_approx_test.cpp`: tests the AAA utility.

Representative full regression inputs outside this directory are in `tests/input_files/` and use the `Reduced_Lung` problem type. Do not confuse those with `Reduced_Lung_1D_Pipe_Flow` inputs, which belong to the ignored `1d_pipe_flow` path.

This documentation was created without running the test suite. For code changes in this directory, run the reduced-lung unit/regression tests appropriate to the local 4C build setup.

## Assumptions and Conventions

Important current assumptions:

- Input topology node ids and element ids are 1-based in input fields and converted to 0-based internal ids.
- Element node ordering is directional: `element_nodes[0]` is inlet/parent and `element_nodes[1]` is outlet/child.
- Junction creation only considers the outlet node of locally owned elements.
- Internal junctions support only connections with two adjacent elements and bifurcations with three adjacent elements. More than three adjacent elements throws.
- Boundary conditions can only be applied to nodes adjacent to exactly one element.
- There can be at most one boundary condition per `(node, type)`.
- Pressure dofs are always offsets `0` and `1` within an element.
- Inlet flow dof is offset `2`; outlet flow dof is the last dof of the element.
- Rigid airway elements have three dofs and one state equation.
- Kelvin-Voigt airway elements have four dofs and two state equations.
- Terminal-unit elements have three dofs and one state equation.
- Model containers group elements by model type combination, so per-element arrays inside a model block must stay aligned.
- Several Jacobian assemblers rely on persistent sparse-matrix state: first call inserts sparsity, later calls replace nonlinear entries or skip constant blocks.
- The current implementation uses `Core::LinAlg::export_to()` to synchronize NOX trial vectors to owned and locally relevant dof vectors.
- Runtime output vectors are element fields and are initialized to `NaN` for unsupported fields/elements.

Dependencies used directly by this directory include 4C core FE/discretization, MPI communication helpers, rebalance utilities, linear algebra maps/vectors/sparse matrices, the global problem/parameter system, the NOX nonlinear solver adapter, Teuchos parameter lists, the function manager/time functions for boundary conditions, and visualization output writers.

## Warnings for Future Agents

- This is baseline documentation for the current pre-thesis implementation. Do not treat the NOX/generic solver workflow as the intended final design.
- Do not implement the future tree-based solver while updating this baseline documentation.
- Keep `src/1d_pipe_flow/` separate unless a task explicitly asks about the 1D pipe-flow module.
- Be careful when changing dof ordering, row ordering, or equation counts. Many helper functions assume fixed offsets and ordering.
- Be careful with the current `x_`/`row_map_` versus `dofs_`/`locally_owned_dof_map_` relationship in `ReducedLungSimulation::build_linear_system_and_solver()` and `NoxSolver::sync_state_from_x()`.
- Be careful with matrix fill state. Some Jacobian assemblers insert only before `SparseMatrix::complete()` and later replace entries, while junction and boundary Jacobians intentionally skip work after the matrix is filled.
- Be careful with node ordering in input topology. Wrong directionality can create wrong parent/child interpretation or trigger duplicate-junction assertions.
- Be careful with boundary condition ownership in parallel runs. Boundary conditions are owned only by the rank that owns the attached element.
- Be careful when adding model variants. Registries group elements by variant pair, and every new model must provide parameter appending, residual/Jacobian callbacks, state updates, output behavior, tests, and compatible dof/state counts.
- The input field API uses zero-based access after conversion even when input maps are documented as 1-based. Existing tests are the best reference for expected behavior.
- `restart_every` is parsed but not used by the current `ReducedLungSimulation` time loop.
- Some test helper discretizations use legacy `Discret::Elements::RedAirway`, while the production topology builder uses lightweight geometry-only line elements. Do not assume the test helper element type is the production element type.
