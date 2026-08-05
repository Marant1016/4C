// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_rebalance.hpp"
#include "4C_reduced_lung_airways.hpp"
#include "4C_reduced_lung_boundary_conditions.hpp"
#include "4C_reduced_lung_helpers.hpp"
#include "4C_reduced_lung_junctions.hpp"
#include "4C_reduced_lung_linear_solver.hpp"
#include "4C_reduced_lung_newton_solver.hpp"
#include "4C_reduced_lung_solver_profile.hpp"
#include "4C_reduced_lung_terminal_unit.hpp"
#include "4C_reduced_lung_tree_linear_solver.hpp"
#include "4C_reduced_lung_tree_linearization.hpp"
#include "4C_reduced_lung_tree_metadata.hpp"
#include "4C_utils_function_manager.hpp"
#include "4C_utils_function_of_time.hpp"

#include <benchmark/benchmark.h>
#include <mpi.h>
#include <Teuchos_ParameterList.hpp>

#include <algorithm>
#include <any>
#include <chrono>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{
  using namespace FourC;
  using namespace FourC::ReducedLung;

  using Clock = std::chrono::steady_clock;
  using ElementType = ReducedLungParameters::LungTree::ElementType;
  using BoundaryType = ReducedLungParameters::BoundaryConditions::Type;
  using ResistanceType = ReducedLungParameters::LungTree::Airways::FlowModel::ResistanceType;
  using WallModelType = ReducedLungParameters::LungTree::Airways::WallModelType;
  using RheologyType =
      ReducedLungParameters::LungTree::TerminalUnits::RheologicalModel::RheologicalModelType;
  using ElasticityType =
      ReducedLungParameters::LungTree::TerminalUnits::ElasticityModel::ElasticityModelType;

  enum class SolverKind
  {
    Nox,
    NewtonSparse,
    NewtonTree,
  };

  enum class BenchmarkCase
  {
    SingleTerminalUnit,
    SerialAirways,
    BalancedAirways,
  };

  double elapsed_seconds(const Clock::time_point start)
  {
    return std::chrono::duration<double>(Clock::now() - start).count();
  }

  bool ensure_serial_benchmark(benchmark::State& state)
  {
    int comm_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    if (comm_size != 1)
    {
      state.SkipWithError("reduced_lung tree benchmarks are serial-only");
      return false;
    }
    return true;
  }

  Core::Utils::FunctionManager make_function_manager(
      const std::vector<std::string>& function_definitions)
  {
    Core::Utils::FunctionManager function_manager;
    std::vector<std::any> functions;
    functions.reserve(function_definitions.size());
    for (const auto& definition : function_definitions)
    {
      functions.emplace_back(std::shared_ptr<Core::Utils::FunctionOfTime>(
          std::make_shared<Core::Utils::SymbolicFunctionOfTime>(
              std::vector<std::string>{definition},
              std::vector<std::shared_ptr<Core::Utils::FunctionVariable>>{})));
    }
    function_manager.set_functions(functions);
    return function_manager;
  }

  ReducedLungParameters::Dynamics make_dynamics(double dt)
  {
    return ReducedLungParameters::Dynamics{
        .time_increment = dt,
        .number_of_steps = 1,
        .restart_every = 1,
        .results_every = 1,
        .linear_solver = 1,
        .max_nonlinear_iterations = 10,
        .nonlinear_residual_tolerance = 1.0e-8,
        .nonlinear_increment_tolerance = 1.0e-10,
    };
  }

  void set_common_air_properties(ReducedLungParameters& params)
  {
    params.air_properties = {
        .density = 1.176e-06,
        .dynamic_viscosity = 1.79105e-05,
    };
  }

  void set_linear_rigid_airway_model(
      ReducedLungParameters& params, const std::unordered_map<int, double>& radii)
  {
    params.lung_tree.airways.radius = Core::IO::InputField<double>(radii);
    params.lung_tree.airways.flow_model.resistance_type =
        Core::IO::InputField<ResistanceType>(ResistanceType::Linear);
    params.lung_tree.airways.flow_model.include_inertia = Core::IO::InputField<bool>(false);
    params.lung_tree.airways.wall_model_type =
        Core::IO::InputField<WallModelType>(WallModelType::Rigid);
  }

  void set_linear_terminal_unit_model(ReducedLungParameters& params)
  {
    params.lung_tree.terminal_units.rheological_model.rheological_model_type =
        Core::IO::InputField<RheologyType>(RheologyType::KelvinVoigt);
    params.lung_tree.terminal_units.rheological_model.kelvin_voigt.viscosity_kelvin_voigt_eta =
        Core::IO::InputField<double>(1.0);
    params.lung_tree.terminal_units.elasticity_model.elasticity_model_type =
        Core::IO::InputField<ElasticityType>(ElasticityType::Linear);
    params.lung_tree.terminal_units.elasticity_model.linear.elasticity_e =
        Core::IO::InputField<double>(1.0);
  }

  ReducedLungParameters make_single_terminal_unit_parameters(double dt)
  {
    ReducedLungParameters params{};
    set_common_air_properties(params);
    params.dynamics = make_dynamics(dt);

    params.lung_tree.topology.num_nodes = 2;
    params.lung_tree.topology.num_elements = 1;
    params.lung_tree.topology.node_coordinates =
        Core::IO::InputField<std::vector<double>>(std::unordered_map<int, std::vector<double>>{
            {1, {0.0, 0.0, 0.0}},
            {2, {1.0, 0.0, 0.0}},
        });
    params.lung_tree.topology.element_nodes =
        Core::IO::InputField<std::vector<int>>(std::unordered_map<int, std::vector<int>>{
            {1, {1, 2}},
        });
    params.lung_tree.element_type = Core::IO::InputField<ElementType>(ElementType::TerminalUnit);
    params.lung_tree.generation = Core::IO::InputField<int>(-1);
    set_linear_terminal_unit_model(params);

    params.boundary_conditions.num_conditions = 2;
    params.boundary_conditions.bc_type = Core::IO::InputField<BoundaryType>(BoundaryType::Pressure);
    params.boundary_conditions.node_id =
        Core::IO::InputField<int>(std::unordered_map<int, int>{{1, 1}, {2, 2}});
    params.boundary_conditions.value_source =
        ReducedLungParameters::BoundaryConditions::ValueSource::bc_function_id;
    params.boundary_conditions.function_id =
        Core::IO::InputField<int>(std::unordered_map<int, int>{{1, 1}, {2, 2}});

    return params;
  }

  ReducedLungParameters make_serial_airway_parameters(int element_count, double dt)
  {
    ReducedLungParameters params{};
    set_common_air_properties(params);
    params.dynamics = make_dynamics(dt);

    std::unordered_map<int, std::vector<double>> coordinates;
    std::unordered_map<int, std::vector<int>> element_nodes;
    std::unordered_map<int, int> generation;
    std::unordered_map<int, double> radii;
    for (int node = 1; node <= element_count + 1; ++node)
    {
      coordinates[node] = {static_cast<double>(node - 1), 0.0, 0.0};
    }
    for (int element = 1; element <= element_count; ++element)
    {
      element_nodes[element] = {element, element + 1};
      generation[element] = element - 1;
      radii[element] = 1.0 / (1.0 + 0.03 * static_cast<double>(element - 1));
    }

    params.lung_tree.topology.num_nodes = element_count + 1;
    params.lung_tree.topology.num_elements = element_count;
    params.lung_tree.topology.node_coordinates =
        Core::IO::InputField<std::vector<double>>(coordinates);
    params.lung_tree.topology.element_nodes = Core::IO::InputField<std::vector<int>>(element_nodes);
    params.lung_tree.element_type = Core::IO::InputField<ElementType>(ElementType::Airway);
    params.lung_tree.generation = Core::IO::InputField<int>(generation);
    set_linear_rigid_airway_model(params, radii);
    set_linear_terminal_unit_model(params);

    params.boundary_conditions.num_conditions = 2;
    params.boundary_conditions.bc_type = Core::IO::InputField<BoundaryType>(BoundaryType::Pressure);
    params.boundary_conditions.node_id =
        Core::IO::InputField<int>(std::unordered_map<int, int>{{1, 1}, {2, element_count + 1}});
    params.boundary_conditions.value_source =
        ReducedLungParameters::BoundaryConditions::ValueSource::bc_function_id;
    params.boundary_conditions.function_id =
        Core::IO::InputField<int>(std::unordered_map<int, int>{{1, 1}, {2, 2}});

    return params;
  }

  ReducedLungParameters make_balanced_airway_parameters(int levels, double dt)
  {
    ReducedLungParameters params{};
    set_common_air_properties(params);
    params.dynamics = make_dynamics(dt);

    std::unordered_map<int, std::vector<double>> coordinates{{1, {0.0, 0.0, 0.0}}};
    std::unordered_map<int, std::vector<int>> element_nodes;
    std::unordered_map<int, int> generation;
    std::unordered_map<int, double> radii;
    std::vector<int> leaf_nodes;
    int next_node = 1;
    int next_element = 1;

    std::function<void(int, int, double)> add_subtree = [&](int inlet_node, int depth, double y)
    {
      const int element = next_element++;
      const int outlet_node = ++next_node;
      coordinates[outlet_node] = {static_cast<double>(depth + 1), y, 0.0};
      element_nodes[element] = {inlet_node, outlet_node};
      generation[element] = depth;
      radii[element] = 1.0 / (1.0 + 0.04 * static_cast<double>(depth));

      if (depth + 1 == levels)
      {
        leaf_nodes.push_back(outlet_node);
        return;
      }

      const double offset = std::pow(0.5, static_cast<double>(depth));
      add_subtree(outlet_node, depth + 1, y + offset);
      add_subtree(outlet_node, depth + 1, y - offset);
    };
    add_subtree(1, 0, 0.0);

    params.lung_tree.topology.num_nodes = next_node;
    params.lung_tree.topology.num_elements = next_element - 1;
    params.lung_tree.topology.node_coordinates =
        Core::IO::InputField<std::vector<double>>(coordinates);
    params.lung_tree.topology.element_nodes = Core::IO::InputField<std::vector<int>>(element_nodes);
    params.lung_tree.element_type = Core::IO::InputField<ElementType>(ElementType::Airway);
    params.lung_tree.generation = Core::IO::InputField<int>(generation);
    set_linear_rigid_airway_model(params, radii);
    set_linear_terminal_unit_model(params);

    std::unordered_map<int, int> boundary_nodes{{1, 1}};
    std::unordered_map<int, int> function_ids{{1, 1}};
    int condition = 2;
    for (const int leaf_node : leaf_nodes)
    {
      boundary_nodes[condition] = leaf_node;
      function_ids[condition] = 2;
      ++condition;
    }
    params.boundary_conditions.num_conditions = static_cast<int>(boundary_nodes.size());
    params.boundary_conditions.bc_type = Core::IO::InputField<BoundaryType>(BoundaryType::Pressure);
    params.boundary_conditions.node_id = Core::IO::InputField<int>(boundary_nodes);
    params.boundary_conditions.value_source =
        ReducedLungParameters::BoundaryConditions::ValueSource::bc_function_id;
    params.boundary_conditions.function_id = Core::IO::InputField<int>(function_ids);

    return params;
  }

  ReducedLungParameters make_benchmark_parameters(BenchmarkCase benchmark_case)
  {
    switch (benchmark_case)
    {
      case BenchmarkCase::SingleTerminalUnit:
        return make_single_terminal_unit_parameters(0.1);
      case BenchmarkCase::SerialAirways:
        return make_serial_airway_parameters(9, 0.1);
      case BenchmarkCase::BalancedAirways:
        return make_balanced_airway_parameters(4, 0.1);
    }
    return make_serial_airway_parameters(3, 0.1);
  }

  std::string benchmark_case_name(BenchmarkCase benchmark_case)
  {
    switch (benchmark_case)
    {
      case BenchmarkCase::SingleTerminalUnit:
        return "single_terminal_unit";
      case BenchmarkCase::SerialAirways:
        return "serial_airways";
      case BenchmarkCase::BalancedAirways:
        return "balanced_airways";
    }
    return "unknown";
  }

  struct BenchmarkFixture
  {
    ReducedLungParameters params;
    Core::Utils::FunctionManager function_manager;
    Core::FE::Discretization discretization;
    Airways::AirwayContainer airways;
    TerminalUnits::TerminalUnitContainer terminal_units;
    std::map<int, int> dof_per_ele;
    int n_airways = 0;
    int n_terminal_units = 0;
    std::map<int, int> first_global_dof_of_ele;
    std::map<int, int> global_dof_per_ele;
    std::map<int, std::vector<int>> global_ele_ids_per_node;
    BoundaryConditions::BoundaryConditionContainer boundary_conditions;
    Junctions::ConnectionData connections;
    Junctions::BifurcationData bifurcations;
    std::unique_ptr<Core::LinAlg::Map> locally_owned_dof_map;
    std::unique_ptr<Core::LinAlg::Map> row_map;
    std::unique_ptr<Core::LinAlg::Map> locally_relevant_dof_map;
    std::unique_ptr<Core::LinAlg::Vector<double>> dofs;
    std::unique_ptr<Core::LinAlg::Vector<double>> locally_relevant_dofs;
    std::unique_ptr<Core::LinAlg::Vector<double>> x;
    std::unique_ptr<Core::LinAlg::SparseMatrix> sysmat;
    ReducedLungAssemblyPipeline assembly_pipeline;
    std::optional<ReducedLungTreeMetadata> tree_metadata;
    Teuchos::ParameterList solver_params;
    std::function<const Teuchos::ParameterList&(int)> solver_params_callback;

    BenchmarkFixture(const std::string& name, ReducedLungParameters input_params)
        : params(std::move(input_params)),
          function_manager(make_function_manager({"0.5 * t", "0.0"})),
          discretization(name, MPI_COMM_WORLD, 3)
    {
      initialize();
    }

    void initialize()
    {
      Core::Rebalance::RebalanceParameters rebalance_parameters;
      build_discretization_from_topology(
          discretization, params.lung_tree.topology, rebalance_parameters);
      discretization.fill_complete();

      create_local_element_models(discretization, params, airways, terminal_units, dof_per_ele,
          n_airways, n_terminal_units);
      create_global_dof_maps(
          dof_per_ele, MPI_COMM_WORLD, global_dof_per_ele, first_global_dof_of_ele);
      assign_global_dof_ids_to_models(first_global_dof_of_ele, airways, terminal_units);
      TerminalUnits::create_evaluators(terminal_units);
      Airways::create_evaluators(airways);

      global_ele_ids_per_node = create_global_ele_ids_per_node(discretization, MPI_COMM_WORLD);
      BoundaryConditions::create_boundary_conditions(discretization, params,
          global_ele_ids_per_node, global_dof_per_ele, first_global_dof_of_ele, function_manager,
          boundary_conditions);
      BoundaryConditions::create_evaluators(boundary_conditions);
      Junctions::create_junctions(discretization, global_ele_ids_per_node, global_dof_per_ele,
          first_global_dof_of_ele, connections, bifurcations);

      int n_local_equations = 0;
      Airways::assign_local_equation_ids(airways, n_local_equations);
      TerminalUnits::assign_local_equation_ids(terminal_units, n_local_equations);
      Junctions::assign_junction_local_equation_ids(connections, bifurcations, n_local_equations);
      BoundaryConditions::assign_local_equation_ids(boundary_conditions, n_local_equations);

      locally_owned_dof_map = std::make_unique<Core::LinAlg::Map>(
          create_domain_map(MPI_COMM_WORLD, airways, terminal_units));
      row_map = std::make_unique<Core::LinAlg::Map>(create_row_map(
          MPI_COMM_WORLD, airways, terminal_units, connections, bifurcations, boundary_conditions));
      locally_relevant_dof_map = std::make_unique<Core::LinAlg::Map>(
          create_column_map(MPI_COMM_WORLD, airways, terminal_units, global_dof_per_ele,
              first_global_dof_of_ele, connections, bifurcations, boundary_conditions));

      Junctions::assign_junction_global_equation_ids(*row_map, connections, bifurcations);
      BoundaryConditions::assign_global_equation_ids(*row_map, boundary_conditions);
      Airways::assign_local_dof_ids(*locally_relevant_dof_map, airways);
      TerminalUnits::assign_local_dof_ids(*locally_relevant_dof_map, terminal_units);
      Junctions::assign_junction_local_dof_ids(
          *locally_relevant_dof_map, connections, bifurcations);
      BoundaryConditions::assign_local_dof_ids(*locally_relevant_dof_map, boundary_conditions);

      dofs = std::make_unique<Core::LinAlg::Vector<double>>(*locally_owned_dof_map, true);
      locally_relevant_dofs =
          std::make_unique<Core::LinAlg::Vector<double>>(*locally_relevant_dof_map, true);
      x = std::make_unique<Core::LinAlg::Vector<double>>(*row_map, true);
      sysmat = std::make_unique<Core::LinAlg::SparseMatrix>(*row_map, *locally_relevant_dof_map, 4);

      solver_params.set("SOLVER", Core::LinearSolver::SolverType::UMFPACK);
      solver_params.set("NAME", "Reduced_lung_benchmark_solver");
      solver_params_callback = [this](int) -> const Teuchos::ParameterList&
      { return solver_params; };
      assembly_pipeline = create_default_reduced_lung_assembly_pipeline(
          airways, terminal_units, connections, bifurcations, boundary_conditions);
    }

    void seed_state(double scale)
    {
      for (const auto& airway_model : airways.models)
      {
        const auto& data = airway_model.data;
        for (std::size_t i = 0; i < data.number_of_elements(); ++i)
        {
          x->replace_global_value(data.gid_p1[i], scale * (1.0 + 0.02 * static_cast<double>(i)));
          x->replace_global_value(data.gid_p2[i], scale * (0.5 + 0.01 * static_cast<double>(i)));
          x->replace_global_value(data.gid_q1[i], scale * (40.0 + 2.0 * static_cast<double>(i)));
          if (i < data.gid_q2.size())
          {
            x->replace_global_value(data.gid_q2[i], scale * (30.0 + 1.5 * static_cast<double>(i)));
          }
        }
      }

      for (const auto& terminal_unit_model : terminal_units.models)
      {
        const auto& data = terminal_unit_model.data;
        for (std::size_t i = 0; i < data.number_of_elements(); ++i)
        {
          x->replace_global_value(data.gid_p1[i], scale * (0.8 + 0.1 * static_cast<double>(i)));
          x->replace_global_value(data.gid_p2[i], scale * (0.2 + 0.05 * static_cast<double>(i)));
          x->replace_global_value(data.gid_q[i], scale * (0.05 + 0.01 * static_cast<double>(i)));
        }
      }
    }

    void sync_state_from_x()
    {
      Core::LinAlg::export_to(*x, *dofs);
      Core::LinAlg::export_to(*dofs, *locally_relevant_dofs);
      for (const auto& update_state : assembly_pipeline.state_updaters)
      {
        update_state(*locally_relevant_dofs, params.dynamics.time_increment);
      }
    }

    Core::LinAlg::Vector<double> assemble_residual(double current_time)
    {
      Core::LinAlg::Vector<double> residual(*row_map, true);
      residual.put_scalar(0.0);
      for (const auto& assemble_residual_callback : assembly_pipeline.residual_assemblers)
      {
        assemble_residual_callback(
            residual, *locally_relevant_dofs, current_time, params.dynamics.time_increment);
      }
      return residual;
    }

    void assemble_jacobian(double current_time)
    {
      for (const auto& assemble_jacobian_callback : assembly_pipeline.jacobian_assemblers)
      {
        assemble_jacobian_callback(
            *sysmat, *locally_relevant_dofs, current_time, params.dynamics.time_increment);
      }
      if (!sysmat->filled()) sysmat->complete();
    }

    TreeLinearization assemble_tree_linearization(double current_time)
    {
      TreeLinearization linearization(
          row_map->num_my_elements(), locally_relevant_dof_map->num_my_elements());
      for (const auto& tree_linearization_assembler :
          assembly_pipeline.tree_linearization_assemblers)
      {
        tree_linearization_assembler.callback(
            linearization, *locally_relevant_dofs, current_time, params.dynamics.time_increment);
      }
      return linearization;
    }

    ReducedLungTreeMetadata build_tree_metadata() const
    {
      return build_reduced_lung_tree_metadata(ReducedLungTreeMetadataContext{
          .parameters = params,
          .first_global_dof_of_ele = first_global_dof_of_ele,
          .global_dof_per_ele = global_dof_per_ele,
          .airways = airways,
          .terminal_units = terminal_units,
          .connections = connections,
          .bifurcations = bifurcations,
          .boundary_conditions = boundary_conditions,
          .row_map = *row_map,
          .locally_relevant_dof_map = *locally_relevant_dof_map,
      });
    }

    std::unique_ptr<NoxSolver> create_nox_solver(NoxSolverProfile* profile)
    {
      const NoxSolverContext context{
          .comm = MPI_COMM_WORLD,
          .dynamics = params.dynamics,
          .linear_solver_parameters = solver_params,
          .solver_params_callback = solver_params_callback,
          .assembly_pipeline = assembly_pipeline,
          .dofs = *dofs,
          .locally_relevant_dofs = *locally_relevant_dofs,
          .x = *x,
          .jacobian = *sysmat,
          .profile = profile,
      };
      return std::make_unique<NoxSolver>(context);
    }

    std::unique_ptr<NewtonSolver> create_newton_solver(
        const std::shared_ptr<NewtonLinearSolver>& linear_solver, NewtonSolverProfile* profile)
    {
      const NewtonSolverContext context{
          .dynamics = params.dynamics,
          .linear_solver = linear_solver,
          .assembly_pipeline = assembly_pipeline,
          .dofs = *dofs,
          .locally_relevant_dofs = *locally_relevant_dofs,
          .x = *x,
          .jacobian = *sysmat,
          .profile = profile,
      };
      return std::make_unique<NewtonSolver>(context);
    }

    std::unique_ptr<NewtonSolver> create_sparse_newton_solver(
        NewtonSolverProfile* newton_profile, SparseNewtonLinearSolverProfile* sparse_profile)
    {
      auto linear_solver =
          std::make_shared<SparseNewtonLinearSolver>(SparseNewtonLinearSolverContext{
              .comm = MPI_COMM_WORLD,
              .linear_solver_parameters = solver_params,
              .solver_params_callback = solver_params_callback,
              .correction_map = *row_map,
              .profile = sparse_profile,
          });
      return create_newton_solver(linear_solver, newton_profile);
    }

    std::unique_ptr<NewtonSolver> create_tree_newton_solver(
        NewtonSolverProfile* newton_profile, TreeNewtonLinearSolverProfile* tree_profile)
    {
      tree_metadata = build_tree_metadata();
      auto linear_solver = std::make_shared<TreeNewtonLinearSolver>(TreeNewtonLinearSolverContext{
          .tree_metadata = *tree_metadata,
          .pivot_tolerance = 1.0e-12,
          .coefficient_source = TreeNewtonLinearSolverCoefficientSource::StructuredTreeBlocks,
          .profile = tree_profile,
      });
      return create_newton_solver(linear_solver, newton_profile);
    }

    double residual_norm(double current_time)
    {
      auto residual = assemble_residual(current_time);
      double norm = 0.0;
      residual.norm_2(&norm);
      return norm;
    }
  };

  void set_common_counters(benchmark::State& state, const BenchmarkFixture& fixture)
  {
    state.counters["elements"] =
        static_cast<double>(fixture.params.lung_tree.topology.num_elements);
    state.counters["dofs"] = static_cast<double>(fixture.locally_owned_dof_map->num_my_elements());
    state.counters["equations"] = static_cast<double>(fixture.row_map->num_my_elements());
  }

  void set_newton_profile_counters(
      benchmark::State& state, const NewtonSolverProfile& profile, double iterations)
  {
    state.counters["newton_total_s"] = profile.total_solve_time / iterations;
    state.counters["state_sync_s"] = profile.state_sync_time / iterations;
    state.counters["residual_s"] = profile.residual_assembly_time / iterations;
    state.counters["sparse_assembly_s"] = profile.sparse_jacobian_assembly_time / iterations;
    state.counters["sparse_complete_s"] = profile.sparse_jacobian_complete_time / iterations;
    state.counters["tree_assembly_s"] =
        profile.structured_tree_linearization_assembly_time / iterations;
    state.counters["tree_assembly_clear_s"] = profile.tree_linearization_clear_time / iterations;
    state.counters["tree_assembly_airways_s"] = profile.tree_linearization_airway_time / iterations;
    state.counters["tree_assembly_terminal_units_s"] =
        profile.tree_linearization_terminal_unit_time / iterations;
    state.counters["tree_assembly_junctions_s"] =
        profile.tree_linearization_junction_time / iterations;
    state.counters["tree_assembly_boundary_conditions_s"] =
        profile.tree_linearization_boundary_condition_time / iterations;
    state.counters["tree_assembly_other_s"] = profile.tree_linearization_other_time / iterations;
    state.counters["tree_assembly_solver_update_s"] =
        profile.tree_linearization_solver_update_time / iterations;
    state.counters["linear_solve_s"] = profile.linear_solve_time / iterations;
  }

  void set_nox_profile_counters(
      benchmark::State& state, const NoxSolverProfile& profile, double iterations)
  {
    state.counters["nox_total_s"] = profile.total_solve_time / iterations;
    state.counters["state_sync_s"] = profile.state_sync_time / iterations;
    state.counters["residual_s"] = profile.residual_assembly_time / iterations;
    state.counters["sparse_assembly_s"] = profile.sparse_jacobian_assembly_time / iterations;
    state.counters["sparse_complete_s"] = profile.sparse_jacobian_complete_time / iterations;
    state.counters["residual_evals"] =
        static_cast<double>(profile.residual_evaluation_count) / iterations;
    state.counters["jacobian_evals"] =
        static_cast<double>(profile.jacobian_evaluation_count) / iterations;
  }

  void set_tree_profile_counters(
      benchmark::State& state, const TreeNewtonLinearSolverProfile& profile, double iterations)
  {
    state.counters["tree_solve_s"] = profile.total_solve_time / iterations;
    state.counters["tree_bottom_up_s"] = profile.bottom_up_time / iterations;
    state.counters["tree_top_down_s"] = profile.top_down_time / iterations;
    state.counters["tree_dense_s"] = profile.dense_solve_time / iterations;
    state.counters["tree_lookup_s"] = profile.coefficient_lookup_time / iterations;
    state.counters["tree_lookups"] =
        static_cast<double>(profile.coefficient_lookup_count) / iterations;
    state.counters["tree_dense_solves"] =
        static_cast<double>(profile.dense_solve_count) / iterations;
    state.counters["tree_simd_groups"] = static_cast<double>(profile.simd_group_count) / iterations;
    state.counters["tree_simd_lanes"] = static_cast<double>(profile.simd_lane_count) / iterations;
    state.counters["tree_scalar_groups"] =
        static_cast<double>(profile.scalar_group_count) / iterations;
    state.counters["tree_scalar_tail_lanes"] =
        static_cast<double>(profile.scalar_tail_lane_count) / iterations;
    state.counters["tree_dense_fallbacks"] =
        static_cast<double>(profile.dense_fallback_count) / iterations;
    state.counters["tree_unsupported_fallbacks"] =
        static_cast<double>(profile.unsupported_block_fallback_count) / iterations;
    state.counters["tree_workspace_dofs"] = static_cast<double>(profile.total_local_block_dofs);
    state.counters["tree_max_block"] = static_cast<double>(profile.max_local_block_size);
  }

  template <SolverKind solver_kind, BenchmarkCase benchmark_case>
  static void full_nonlinear_solve(benchmark::State& state)
  {
    if (!ensure_serial_benchmark(state)) return;

    NoxSolverProfile nox_profile;
    NewtonSolverProfile newton_profile;
    SparseNewtonLinearSolverProfile sparse_profile;
    TreeNewtonLinearSolverProfile tree_profile;
    std::uint64_t nonlinear_iterations = 0;
    double final_residual_norm = 0.0;
    double element_count = 0.0;
    double dof_count = 0.0;
    double equation_count = 0.0;

    for (auto _ : state)
    {
      state.PauseTiming();
      unsigned int iterations = 0;
      double iteration_residual_norm = 0.0;
      double iteration_element_count = 0.0;
      double iteration_dof_count = 0.0;
      double iteration_equation_count = 0.0;

      {
        auto params = make_benchmark_parameters(benchmark_case);
        const double current_time = params.dynamics.time_increment;
        BenchmarkFixture fixture(benchmark_case_name(benchmark_case), std::move(params));
        fixture.seed_state(0.01);

        if constexpr (solver_kind == SolverKind::Nox)
        {
          auto solver = fixture.create_nox_solver(&nox_profile);
          state.ResumeTiming();
          iterations = solver->solve(current_time);
          benchmark::DoNotOptimize(iterations);
          state.PauseTiming();
        }
        else if constexpr (solver_kind == SolverKind::NewtonSparse)
        {
          auto solver = fixture.create_sparse_newton_solver(&newton_profile, &sparse_profile);
          state.ResumeTiming();
          iterations = solver->solve(current_time);
          benchmark::DoNotOptimize(iterations);
          state.PauseTiming();
        }
        else
        {
          auto solver = fixture.create_tree_newton_solver(&newton_profile, &tree_profile);
          state.ResumeTiming();
          iterations = solver->solve(current_time);
          benchmark::DoNotOptimize(iterations);
          state.PauseTiming();
        }

        iteration_residual_norm = fixture.residual_norm(current_time);
        iteration_element_count =
            static_cast<double>(fixture.params.lung_tree.topology.num_elements);
        iteration_dof_count = static_cast<double>(fixture.locally_owned_dof_map->num_my_elements());
        iteration_equation_count = static_cast<double>(fixture.row_map->num_my_elements());
      }

      nonlinear_iterations += iterations;
      final_residual_norm += iteration_residual_norm;
      element_count = iteration_element_count;
      dof_count = iteration_dof_count;
      equation_count = iteration_equation_count;
      state.ResumeTiming();
    }

    const double benchmark_iterations = static_cast<double>(state.iterations());
    state.counters["nonlinear_iterations"] =
        static_cast<double>(nonlinear_iterations) / benchmark_iterations;
    state.counters["final_residual"] = final_residual_norm / benchmark_iterations;
    state.counters["elements"] = element_count;
    state.counters["dofs"] = dof_count;
    state.counters["equations"] = equation_count;
    if constexpr (solver_kind == SolverKind::Nox)
    {
      set_nox_profile_counters(state, nox_profile, benchmark_iterations);
    }
    else
    {
      set_newton_profile_counters(state, newton_profile, benchmark_iterations);
      if constexpr (solver_kind == SolverKind::NewtonSparse)
      {
        state.counters["sparse_solve_s"] = sparse_profile.solve_time / benchmark_iterations;
      }
      else
      {
        set_tree_profile_counters(state, tree_profile, benchmark_iterations);
      }
    }
  }

  static void assembly_phases_balanced_airways(benchmark::State& state)
  {
    if (!ensure_serial_benchmark(state)) return;

    BenchmarkFixture fixture(
        "balanced_airway_assembly_phases", make_balanced_airway_parameters(4, 0.1));
    fixture.seed_state(0.01);
    fixture.sync_state_from_x();
    const double current_time = fixture.params.dynamics.time_increment;
    double residual_time = 0.0;
    double sparse_assembly_time = 0.0;
    double sparse_complete_time = 0.0;
    double tree_assembly_time = 0.0;

    for (auto _ : state)
    {
      Core::LinAlg::Vector<double> residual(*fixture.row_map, true);
      residual.put_scalar(0.0);
      auto phase_start = Clock::now();
      for (const auto& assemble_residual_callback : fixture.assembly_pipeline.residual_assemblers)
      {
        assemble_residual_callback(residual, *fixture.locally_relevant_dofs, current_time,
            fixture.params.dynamics.time_increment);
      }
      residual_time += elapsed_seconds(phase_start);
      benchmark::DoNotOptimize(residual.local_length());

      fixture.sysmat = std::make_unique<Core::LinAlg::SparseMatrix>(
          *fixture.row_map, *fixture.locally_relevant_dof_map, 4);
      phase_start = Clock::now();
      for (const auto& assemble_jacobian_callback : fixture.assembly_pipeline.jacobian_assemblers)
      {
        assemble_jacobian_callback(*fixture.sysmat, *fixture.locally_relevant_dofs, current_time,
            fixture.params.dynamics.time_increment);
      }
      sparse_assembly_time += elapsed_seconds(phase_start);

      phase_start = Clock::now();
      fixture.sysmat->complete();
      sparse_complete_time += elapsed_seconds(phase_start);
      benchmark::DoNotOptimize(fixture.sysmat->filled());

      TreeLinearization tree_linearization(
          fixture.row_map->num_my_elements(), fixture.locally_relevant_dof_map->num_my_elements());
      phase_start = Clock::now();
      for (const auto& tree_linearization_assembler :
          fixture.assembly_pipeline.tree_linearization_assemblers)
      {
        tree_linearization_assembler.callback(tree_linearization, *fixture.locally_relevant_dofs,
            current_time, fixture.params.dynamics.time_increment);
      }
      tree_assembly_time += elapsed_seconds(phase_start);
      benchmark::DoNotOptimize(tree_linearization.num_rows());
    }

    const double iterations = static_cast<double>(state.iterations());
    set_common_counters(state, fixture);
    state.counters["residual_s"] = residual_time / iterations;
    state.counters["sparse_assembly_s"] = sparse_assembly_time / iterations;
    state.counters["sparse_complete_s"] = sparse_complete_time / iterations;
    state.counters["tree_assembly_s"] = tree_assembly_time / iterations;
  }

  static void sparse_linear_solve_balanced_airways(benchmark::State& state)
  {
    if (!ensure_serial_benchmark(state)) return;

    const int levels = static_cast<int>(state.range(0));
    BenchmarkFixture fixture(
        "balanced_airway_sparse_linear_solve", make_balanced_airway_parameters(levels, 0.1));
    fixture.seed_state(0.01);
    fixture.sync_state_from_x();
    const double current_time = fixture.params.dynamics.time_increment;
    auto residual = fixture.assemble_residual(current_time);
    fixture.assemble_jacobian(current_time);

    SparseNewtonLinearSolverProfile sparse_profile;
    SparseNewtonLinearSolver sparse_solver(SparseNewtonLinearSolverContext{
        .comm = MPI_COMM_WORLD,
        .linear_solver_parameters = fixture.solver_params,
        .solver_params_callback = fixture.solver_params_callback,
        .correction_map = *fixture.row_map,
        .profile = &sparse_profile,
    });
    Core::LinAlg::Vector<double> delta(*fixture.row_map, true);
    const NewtonLinearSystemMetadata metadata{.current_time = current_time,
        .time_step_size_dt = fixture.params.dynamics.time_increment,
        .nonlinear_iteration = 0};

    for (auto _ : state)
    {
      sparse_solver.solve(*fixture.sysmat, residual, *fixture.x, metadata, delta);
      benchmark::DoNotOptimize(delta.local_length());
    }

    const double iterations = static_cast<double>(state.iterations());
    set_common_counters(state, fixture);
    state.counters["sparse_solve_s"] = sparse_profile.solve_time / iterations;
  }

  static void structured_tree_linear_solve_balanced_airways(benchmark::State& state)
  {
    if (!ensure_serial_benchmark(state)) return;

    const int levels = static_cast<int>(state.range(0));
    BenchmarkFixture fixture(
        "balanced_airway_tree_linear_solve", make_balanced_airway_parameters(levels, 0.1));
    fixture.seed_state(0.01);
    fixture.sync_state_from_x();
    const double current_time = fixture.params.dynamics.time_increment;
    auto residual = fixture.assemble_residual(current_time);
    fixture.assemble_jacobian(current_time);
    auto tree_linearization = fixture.assemble_tree_linearization(current_time);
    const auto tree_metadata = fixture.build_tree_metadata();

    TreeNewtonLinearSolverProfile tree_profile;
    TreeNewtonLinearSolver tree_solver(TreeNewtonLinearSolverContext{
        .tree_metadata = tree_metadata,
        .pivot_tolerance = 1.0e-12,
        .coefficient_source = TreeNewtonLinearSolverCoefficientSource::StructuredTreeBlocks,
        .profile = &tree_profile,
    });
    tree_solver.set_tree_linearization(tree_linearization);
    Core::LinAlg::Vector<double> delta(*fixture.row_map, true);
    const NewtonLinearSystemMetadata metadata{.current_time = current_time,
        .time_step_size_dt = fixture.params.dynamics.time_increment,
        .nonlinear_iteration = 0};

    for (auto _ : state)
    {
      tree_solver.solve(*fixture.sysmat, residual, *fixture.x, metadata, delta);
      benchmark::DoNotOptimize(delta.local_length());
    }

    const double iterations = static_cast<double>(state.iterations());
    set_common_counters(state, fixture);
    set_tree_profile_counters(state, tree_profile, iterations);
  }
}  // namespace

BENCHMARK(full_nonlinear_solve<SolverKind::Nox, BenchmarkCase::SingleTerminalUnit>)
    ->Name("ReducedLung/FullSolve/SingleTerminalUnit/Nox")
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(full_nonlinear_solve<SolverKind::NewtonSparse, BenchmarkCase::SingleTerminalUnit>)
    ->Name("ReducedLung/FullSolve/SingleTerminalUnit/NewtonSparse")
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(full_nonlinear_solve<SolverKind::NewtonTree, BenchmarkCase::SingleTerminalUnit>)
    ->Name("ReducedLung/FullSolve/SingleTerminalUnit/NewtonTree")
    ->Unit(benchmark::kMicrosecond);

BENCHMARK(full_nonlinear_solve<SolverKind::Nox, BenchmarkCase::SerialAirways>)
    ->Name("ReducedLung/FullSolve/SerialAirways/Nox")
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(full_nonlinear_solve<SolverKind::NewtonSparse, BenchmarkCase::SerialAirways>)
    ->Name("ReducedLung/FullSolve/SerialAirways/NewtonSparse")
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(full_nonlinear_solve<SolverKind::NewtonTree, BenchmarkCase::SerialAirways>)
    ->Name("ReducedLung/FullSolve/SerialAirways/NewtonTree")
    ->Unit(benchmark::kMicrosecond);

BENCHMARK(full_nonlinear_solve<SolverKind::Nox, BenchmarkCase::BalancedAirways>)
    ->Name("ReducedLung/FullSolve/BalancedAirways/Nox")
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(full_nonlinear_solve<SolverKind::NewtonSparse, BenchmarkCase::BalancedAirways>)
    ->Name("ReducedLung/FullSolve/BalancedAirways/NewtonSparse")
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(full_nonlinear_solve<SolverKind::NewtonTree, BenchmarkCase::BalancedAirways>)
    ->Name("ReducedLung/FullSolve/BalancedAirways/NewtonTree")
    ->Unit(benchmark::kMicrosecond);

BENCHMARK(assembly_phases_balanced_airways)
    ->Name("ReducedLung/AssemblyPhases/BalancedAirways")
    ->Unit(benchmark::kMicrosecond);

BENCHMARK(sparse_linear_solve_balanced_airways)
    ->Name("ReducedLung/LinearSolve/BalancedAirways/Sparse")
    ->Arg(2)
    ->Arg(3)
    ->Arg(4)
    ->Arg(5)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(structured_tree_linear_solve_balanced_airways)
    ->Name("ReducedLung/LinearSolve/BalancedAirways/StructuredTree")
    ->Arg(2)
    ->Arg(3)
    ->Arg(4)
    ->Arg(5)
    ->Unit(benchmark::kMicrosecond);
