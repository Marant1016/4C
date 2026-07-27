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
#include <cmath>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{
  using namespace FourC;
  using namespace FourC::ReducedLung;

  using BoundaryType = ReducedLungParameters::BoundaryConditions::Type;
  using ElasticityType =
      ReducedLungParameters::LungTree::TerminalUnits::ElasticityModel::ElasticityModelType;
  using ElementType = ReducedLungParameters::LungTree::ElementType;
  using ResistanceType = ReducedLungParameters::LungTree::Airways::FlowModel::ResistanceType;
  using RheologyType =
      ReducedLungParameters::LungTree::TerminalUnits::RheologicalModel::RheologicalModelType;
  using WallModelType = ReducedLungParameters::LungTree::Airways::WallModelType;

  enum class SolverKind
  {
    Nox,
    NewtonSparse,
    NewtonTree,
  };

  int mpi_size(MPI_Comm comm)
  {
    int size = 1;
    MPI_Comm_size(comm, &size);
    return size;
  }

  double mpi_max(double value, MPI_Comm comm)
  {
    double global_value = 0.0;
    MPI_Allreduce(&value, &global_value, 1, MPI_DOUBLE, MPI_MAX, comm);
    return global_value;
  }

  std::uint64_t mpi_sum(std::uint64_t value, MPI_Comm comm)
  {
    auto local_value = static_cast<unsigned long long>(value);
    unsigned long long global_value = 0;
    MPI_Allreduce(&local_value, &global_value, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
    return static_cast<std::uint64_t>(global_value);
  }

  int mpi_max(int value, MPI_Comm comm)
  {
    int global_value = 0;
    MPI_Allreduce(&value, &global_value, 1, MPI_INT, MPI_MAX, comm);
    return global_value;
  }

  unsigned int mpi_max(unsigned int value, MPI_Comm comm)
  {
    unsigned int global_value = 0;
    MPI_Allreduce(&value, &global_value, 1, MPI_UNSIGNED, MPI_MAX, comm);
    return global_value;
  }

  double time_mpi_region(MPI_Comm comm, const std::function<void()>& operation)
  {
    MPI_Barrier(comm);
    const double start = MPI_Wtime();
    operation();
    MPI_Barrier(comm);
    const double local_elapsed = MPI_Wtime() - start;
    return mpi_max(local_elapsed, comm);
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

  void set_linear_terminal_unit_defaults(ReducedLungParameters& params)
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
    set_linear_terminal_unit_defaults(params);

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

  int count_cross_rank_edges(const ReducedLungTreeMetadata& tree_metadata)
  {
    int cross_rank_edges = 0;
    for (const auto& element : tree_metadata.elements)
    {
      if (element.parent_element_index == -1)
      {
        continue;
      }
      const auto& parent =
          tree_metadata.elements[static_cast<std::size_t>(element.parent_element_index)];
      if (parent.owner_rank != element.owner_rank)
      {
        ++cross_rank_edges;
      }
    }
    return cross_rank_edges;
  }

  struct DistributedBenchmarkFixture
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

    DistributedBenchmarkFixture(const std::string& name, ReducedLungParameters input_params)
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
      solver_params.set("NAME", "Reduced_lung_distributed_benchmark_solver");
      solver_params_callback = [this](int) -> const Teuchos::ParameterList&
      { return solver_params; };
      assembly_pipeline = create_default_reduced_lung_assembly_pipeline(
          airways, terminal_units, connections, bifurcations, boundary_conditions);
    }

    void seed_state(double scale)
    {
      x->put_scalar(0.0);
      const auto& solution_map = x->get_map();
      const int global_dofs = locally_owned_dof_map->num_global_elements();
      for (int lid = 0; lid < solution_map.num_my_elements(); ++lid)
      {
        const int gid = solution_map.gid(lid);
        if (gid >= global_dofs)
        {
          continue;
        }

        const double element_index = static_cast<double>(gid / 3);
        double value = scale * (40.0 + 2.0 * element_index);
        if (gid % 3 == 0) value = scale * (1.0 + 0.02 * element_index);
        if (gid % 3 == 1) value = scale * (0.5 + 0.01 * element_index);
        x->replace_global_value(gid, value);
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

    TreeLinearization assemble_tree_linearization(double current_time)
    {
      TreeLinearization linearization(
          row_map->num_my_elements(), locally_relevant_dof_map->num_my_elements());
      for (const auto& assemble_tree_linearization_callback :
          assembly_pipeline.tree_linearization_assemblers)
      {
        assemble_tree_linearization_callback(
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

    std::unique_ptr<NewtonSolver> create_distributed_tree_newton_solver(
        NewtonSolverProfile* newton_profile, TreeNewtonLinearSolverProfile* tree_profile)
    {
      tree_metadata = build_tree_metadata();
      auto linear_solver = std::make_shared<DistributedTreeNewtonLinearSolver>(
          DistributedTreeNewtonLinearSolverContext{
              .tree_metadata = *tree_metadata,
              .locally_relevant_dof_map = *locally_relevant_dof_map,
              .pivot_tolerance = 1.0e-12,
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

  void set_common_counters(benchmark::State& state, const DistributedBenchmarkFixture& fixture,
      const ReducedLungTreeMetadata& tree_metadata)
  {
    state.counters["mpi_ranks"] = static_cast<double>(mpi_size(MPI_COMM_WORLD));
    state.counters["elements"] =
        static_cast<double>(fixture.params.lung_tree.topology.num_elements);
    state.counters["dofs"] =
        static_cast<double>(fixture.locally_owned_dof_map->num_global_elements());
    state.counters["equations"] = static_cast<double>(fixture.row_map->num_global_elements());
    state.counters["cross_rank_edges"] = static_cast<double>(count_cross_rank_edges(tree_metadata));
  }

  void set_newton_profile_counters(
      benchmark::State& state, const NewtonSolverProfile& profile, double iterations)
  {
    state.counters["newton_total_s"] =
        mpi_max(profile.total_solve_time, MPI_COMM_WORLD) / iterations;
    state.counters["state_sync_s"] = mpi_max(profile.state_sync_time, MPI_COMM_WORLD) / iterations;
    state.counters["residual_s"] =
        mpi_max(profile.residual_assembly_time, MPI_COMM_WORLD) / iterations;
    state.counters["sparse_assembly_s"] =
        mpi_max(profile.sparse_jacobian_assembly_time, MPI_COMM_WORLD) / iterations;
    state.counters["sparse_complete_s"] =
        mpi_max(profile.sparse_jacobian_complete_time, MPI_COMM_WORLD) / iterations;
    state.counters["tree_assembly_s"] =
        mpi_max(profile.structured_tree_linearization_assembly_time, MPI_COMM_WORLD) / iterations;
    state.counters["linear_solve_s"] =
        mpi_max(profile.linear_solve_time, MPI_COMM_WORLD) / iterations;
  }

  void set_nox_profile_counters(
      benchmark::State& state, const NoxSolverProfile& profile, double iterations)
  {
    state.counters["nox_total_s"] = mpi_max(profile.total_solve_time, MPI_COMM_WORLD) / iterations;
    state.counters["state_sync_s"] = mpi_max(profile.state_sync_time, MPI_COMM_WORLD) / iterations;
    state.counters["residual_s"] =
        mpi_max(profile.residual_assembly_time, MPI_COMM_WORLD) / iterations;
    state.counters["sparse_assembly_s"] =
        mpi_max(profile.sparse_jacobian_assembly_time, MPI_COMM_WORLD) / iterations;
    state.counters["sparse_complete_s"] =
        mpi_max(profile.sparse_jacobian_complete_time, MPI_COMM_WORLD) / iterations;
    state.counters["residual_evals"] =
        static_cast<double>(
            mpi_max(static_cast<int>(profile.residual_evaluation_count), MPI_COMM_WORLD)) /
        iterations;
    state.counters["jacobian_evals"] =
        static_cast<double>(
            mpi_max(static_cast<int>(profile.jacobian_evaluation_count), MPI_COMM_WORLD)) /
        iterations;
  }

  void set_tree_profile_counters(
      benchmark::State& state, const TreeNewtonLinearSolverProfile& profile, double iterations)
  {
    state.counters["tree_solve_s"] = mpi_max(profile.total_solve_time, MPI_COMM_WORLD) / iterations;
    state.counters["tree_bottom_up_s"] =
        mpi_max(profile.bottom_up_time, MPI_COMM_WORLD) / iterations;
    state.counters["tree_top_down_s"] = mpi_max(profile.top_down_time, MPI_COMM_WORLD) / iterations;
    state.counters["tree_dense_s"] = mpi_max(profile.dense_solve_time, MPI_COMM_WORLD) / iterations;
    state.counters["tree_lookup_s"] =
        mpi_max(profile.coefficient_lookup_time, MPI_COMM_WORLD) / iterations;
    state.counters["tree_comm_s"] =
        mpi_max(profile.communication_time, MPI_COMM_WORLD) / iterations;
    state.counters["tree_dense_solves"] =
        static_cast<double>(mpi_sum(profile.dense_solve_count, MPI_COMM_WORLD)) / iterations;
    state.counters["tree_lookups"] =
        static_cast<double>(mpi_sum(profile.coefficient_lookup_count, MPI_COMM_WORLD)) / iterations;
    state.counters["relation_messages"] =
        static_cast<double>(mpi_sum(profile.boundary_relation_message_count, MPI_COMM_WORLD)) /
        iterations;
    state.counters["pressure_messages"] =
        static_cast<double>(mpi_sum(profile.boundary_pressure_message_count, MPI_COMM_WORLD)) /
        iterations;
    state.counters["scatter_messages"] =
        static_cast<double>(mpi_sum(profile.correction_scatter_message_count, MPI_COMM_WORLD)) /
        iterations;
    state.counters["communication_bytes"] =
        static_cast<double>(mpi_sum(profile.communication_bytes, MPI_COMM_WORLD)) / iterations;
    state.counters["tree_workspace_dofs"] =
        static_cast<double>(mpi_sum(profile.total_local_block_dofs, MPI_COMM_WORLD));
    state.counters["tree_max_block"] =
        static_cast<double>(mpi_max(profile.max_local_block_size, MPI_COMM_WORLD));
  }

  template <SolverKind solver_kind>
  void distributed_full_solve_balanced_airways(benchmark::State& state)
  {
    NoxSolverProfile nox_profile;
    NewtonSolverProfile newton_profile;
    SparseNewtonLinearSolverProfile sparse_profile;
    TreeNewtonLinearSolverProfile tree_profile;
    std::uint64_t nonlinear_iterations = 0;
    double final_residual_norm = 0.0;
    double last_cross_rank_edges = 0.0;
    double last_dofs = 0.0;
    double last_elements = 0.0;
    double last_equations = 0.0;

    for (auto _ : state)
    {
      auto params = make_balanced_airway_parameters(4, 0.1);
      const double current_time = params.dynamics.time_increment;
      DistributedBenchmarkFixture fixture("distributed_balanced_full_solve", std::move(params));
      fixture.seed_state(0.01);
      const auto tree_metadata = fixture.build_tree_metadata();

      unsigned int iterations = 0;
      double elapsed = 0.0;
      if constexpr (solver_kind == SolverKind::Nox)
      {
        auto solver = fixture.create_nox_solver(&nox_profile);
        elapsed =
            time_mpi_region(MPI_COMM_WORLD, [&]() { iterations = solver->solve(current_time); });
      }
      else if constexpr (solver_kind == SolverKind::NewtonSparse)
      {
        auto solver = fixture.create_sparse_newton_solver(&newton_profile, &sparse_profile);
        elapsed =
            time_mpi_region(MPI_COMM_WORLD, [&]() { iterations = solver->solve(current_time); });
      }
      else
      {
        auto solver = fixture.create_distributed_tree_newton_solver(&newton_profile, &tree_profile);
        elapsed =
            time_mpi_region(MPI_COMM_WORLD, [&]() { iterations = solver->solve(current_time); });
      }

      state.SetIterationTime(elapsed);
      nonlinear_iterations += mpi_max(iterations, MPI_COMM_WORLD);
      final_residual_norm += fixture.residual_norm(current_time);
      last_cross_rank_edges = static_cast<double>(count_cross_rank_edges(tree_metadata));
      last_elements = static_cast<double>(fixture.params.lung_tree.topology.num_elements);
      last_dofs = static_cast<double>(fixture.locally_owned_dof_map->num_global_elements());
      last_equations = static_cast<double>(fixture.row_map->num_global_elements());
    }

    const double iterations = static_cast<double>(state.iterations());
    state.counters["mpi_ranks"] = static_cast<double>(mpi_size(MPI_COMM_WORLD));
    state.counters["elements"] = last_elements;
    state.counters["dofs"] = last_dofs;
    state.counters["equations"] = last_equations;
    state.counters["cross_rank_edges"] = last_cross_rank_edges;
    state.counters["nonlinear_iterations"] = static_cast<double>(nonlinear_iterations) / iterations;
    state.counters["final_residual"] = final_residual_norm / iterations;

    if constexpr (solver_kind == SolverKind::Nox)
    {
      set_nox_profile_counters(state, nox_profile, iterations);
    }
    else
    {
      set_newton_profile_counters(state, newton_profile, iterations);
      if constexpr (solver_kind == SolverKind::NewtonSparse)
      {
        state.counters["sparse_solve_s"] =
            mpi_max(sparse_profile.solve_time, MPI_COMM_WORLD) / iterations;
      }
      else
      {
        set_tree_profile_counters(state, tree_profile, iterations);
      }
    }
  }

  void distributed_tree_linear_solve_balanced_airways(benchmark::State& state)
  {
    const int levels = static_cast<int>(state.range(0));
    auto params = make_balanced_airway_parameters(levels, 0.1);
    const double current_time = params.dynamics.time_increment;
    DistributedBenchmarkFixture fixture(
        "distributed_balanced_tree_linear_solve", std::move(params));
    fixture.seed_state(0.01);
    fixture.sync_state_from_x();
    auto residual = fixture.assemble_residual(current_time);
    auto tree_linearization = fixture.assemble_tree_linearization(current_time);
    const auto tree_metadata = fixture.build_tree_metadata();

    TreeNewtonLinearSolverProfile tree_profile;
    DistributedTreeNewtonLinearSolver tree_solver(DistributedTreeNewtonLinearSolverContext{
        .tree_metadata = tree_metadata,
        .locally_relevant_dof_map = *fixture.locally_relevant_dof_map,
        .pivot_tolerance = 1.0e-12,
        .profile = &tree_profile,
    });
    tree_solver.set_tree_linearization(tree_linearization);
    Core::LinAlg::Vector<double> delta(*fixture.row_map, true);
    const NewtonLinearSystemMetadata metadata{.current_time = current_time,
        .time_step_size_dt = fixture.params.dynamics.time_increment,
        .nonlinear_iteration = 0};

    for (auto _ : state)
    {
      const double elapsed = time_mpi_region(MPI_COMM_WORLD,
          [&]() { tree_solver.solve(*fixture.sysmat, residual, *fixture.x, metadata, delta); });
      state.SetIterationTime(elapsed);
    }

    const double iterations = static_cast<double>(state.iterations());
    set_common_counters(state, fixture, tree_metadata);
    set_tree_profile_counters(state, tree_profile, iterations);
  }
}  // namespace

BENCHMARK(distributed_full_solve_balanced_airways<SolverKind::Nox>)
    ->Name("ReducedLung/Distributed/FullSolve/BalancedAirways/Nox")
    ->UseManualTime()
    ->Iterations(10)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(distributed_full_solve_balanced_airways<SolverKind::NewtonSparse>)
    ->Name("ReducedLung/Distributed/FullSolve/BalancedAirways/NewtonSparse")
    ->UseManualTime()
    ->Iterations(10)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(distributed_full_solve_balanced_airways<SolverKind::NewtonTree>)
    ->Name("ReducedLung/Distributed/FullSolve/BalancedAirways/NewtonTree")
    ->UseManualTime()
    ->Iterations(10)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK(distributed_tree_linear_solve_balanced_airways)
    ->Name("ReducedLung/Distributed/LinearSolve/BalancedAirways/StructuredTree")
    ->UseManualTime()
    ->Iterations(20)
    ->Arg(2)
    ->Arg(3)
    ->Arg(4)
    ->Arg(5)
    ->Unit(benchmark::kMicrosecond);
