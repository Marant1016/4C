// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_fem_discretization.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_rebalance.hpp"
#include "4C_reduced_lung_boundary_conditions.hpp"
#include "4C_reduced_lung_helpers.hpp"
#include "4C_reduced_lung_junctions.hpp"
#include "4C_reduced_lung_linear_solver.hpp"
#include "4C_reduced_lung_terminal_unit.hpp"
#include "4C_reduced_lung_tree_linear_solver.hpp"
#include "4C_reduced_lung_tree_linearization.hpp"
#include "4C_reduced_lung_tree_metadata.hpp"
#include "4C_utils_function_manager.hpp"
#include "4C_utils_function_of_time.hpp"

#include <mpi.h>
#include <Teuchos_ParameterList.hpp>

#include <any>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{
  using namespace FourC;
  using namespace FourC::ReducedLung;

  using BoundaryType = ReducedLungParameters::BoundaryConditions::Type;
  using ElementType = ReducedLungParameters::LungTree::ElementType;
  using ElasticityType =
      ReducedLungParameters::LungTree::TerminalUnits::ElasticityModel::ElasticityModelType;
  using ResistanceType = ReducedLungParameters::LungTree::Airways::FlowModel::ResistanceType;
  using RheologyType =
      ReducedLungParameters::LungTree::TerminalUnits::RheologicalModel::RheologicalModelType;
  using WallModelType = ReducedLungParameters::LungTree::Airways::WallModelType;

  Core::Utils::FunctionManager make_function_manager()
  {
    Core::Utils::FunctionManager function_manager;
    std::vector<std::any> functions;
    functions.emplace_back(std::shared_ptr<Core::Utils::FunctionOfTime>(
        std::make_shared<Core::Utils::SymbolicFunctionOfTime>(std::vector<std::string>{"t"},
            std::vector<std::shared_ptr<Core::Utils::FunctionVariable>>{})));
    functions.emplace_back(std::shared_ptr<Core::Utils::FunctionOfTime>(
        std::make_shared<Core::Utils::SymbolicFunctionOfTime>(std::vector<std::string>{"0.0"},
            std::vector<std::shared_ptr<Core::Utils::FunctionVariable>>{})));
    function_manager.set_functions(functions);
    return function_manager;
  }

  ReducedLungParameters make_serial_airway_parameters(double dt)
  {
    ReducedLungParameters params{};
    params.air_properties = {
        .density = 1.176e-06,
        .dynamic_viscosity = 1.79105e-05,
    };
    params.dynamics = ReducedLungParameters::Dynamics{
        .time_increment = dt,
        .number_of_steps = 1,
        .restart_every = 1,
        .results_every = 1,
        .linear_solver = 1,
        .max_nonlinear_iterations = 10,
        .nonlinear_residual_tolerance = 1.0e-8,
        .nonlinear_increment_tolerance = 1.0e-10,
    };

    params.lung_tree.topology.num_nodes = 4;
    params.lung_tree.topology.num_elements = 3;
    params.lung_tree.topology.node_coordinates =
        Core::IO::InputField<std::vector<double>>(std::unordered_map<int, std::vector<double>>{
            {1, {0.0, 0.0, 0.0}},
            {2, {1.0, 0.0, 0.0}},
            {3, {2.0, 0.0, 0.0}},
            {4, {3.0, 0.0, 0.0}},
        });
    params.lung_tree.topology.element_nodes =
        Core::IO::InputField<std::vector<int>>(std::unordered_map<int, std::vector<int>>{
            {1, {1, 2}},
            {2, {2, 3}},
            {3, {3, 4}},
        });
    params.lung_tree.element_type = Core::IO::InputField<ElementType>(ElementType::Airway);
    params.lung_tree.generation =
        Core::IO::InputField<int>(std::unordered_map<int, int>{{1, 0}, {2, 1}, {3, 2}});

    params.lung_tree.airways.radius =
        Core::IO::InputField<double>(std::unordered_map<int, double>{{1, 1.0}, {2, 0.9}, {3, 0.8}});
    params.lung_tree.airways.flow_model.resistance_type =
        Core::IO::InputField<ResistanceType>(ResistanceType::Linear);
    params.lung_tree.airways.flow_model.include_inertia = Core::IO::InputField<bool>(false);
    params.lung_tree.airways.wall_model_type =
        Core::IO::InputField<WallModelType>(WallModelType::Rigid);
    params.lung_tree.terminal_units.rheological_model.rheological_model_type =
        Core::IO::InputField<RheologyType>(RheologyType::KelvinVoigt);
    params.lung_tree.terminal_units.elasticity_model.elasticity_model_type =
        Core::IO::InputField<ElasticityType>(ElasticityType::Linear);
    params.lung_tree.terminal_units.elasticity_model.linear.elasticity_e =
        Core::IO::InputField<double>(1.0);

    params.boundary_conditions.num_conditions = 2;
    params.boundary_conditions.bc_type = Core::IO::InputField<BoundaryType>(BoundaryType::Pressure);
    params.boundary_conditions.node_id =
        Core::IO::InputField<int>(std::unordered_map<int, int>{{1, 1}, {2, 4}});
    params.boundary_conditions.value_source =
        ReducedLungParameters::BoundaryConditions::ValueSource::bc_function_id;
    params.boundary_conditions.function_id =
        Core::IO::InputField<int>(std::unordered_map<int, int>{{1, 1}, {2, 2}});

    return params;
  }

  struct DistributedTreeFixture
  {
    ReducedLungParameters params = make_serial_airway_parameters(0.1);
    Core::Utils::FunctionManager function_manager = make_function_manager();
    Core::FE::Discretization discretization{"distributed_tree_test", MPI_COMM_WORLD, 3};
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
    Teuchos::ParameterList solver_params;
    std::function<const Teuchos::ParameterList&(int)> solver_params_callback;

    DistributedTreeFixture()
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
      solver_params.set("NAME", "Reduced_lung_solver");
      solver_params_callback = [this](int) -> const Teuchos::ParameterList&
      { return solver_params; };
      assembly_pipeline = create_default_reduced_lung_assembly_pipeline(
          airways, terminal_units, connections, bifurcations, boundary_conditions);
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
  };

  void expect_vectors_near(const Core::LinAlg::Vector<double>& expected,
      const Core::LinAlg::Vector<double>& actual, double tolerance)
  {
    ASSERT_EQ(expected.global_length(), actual.global_length());
    Core::LinAlg::Vector<double> actual_on_expected_map(expected.get_map(), true);
    Core::LinAlg::export_to(actual, actual_on_expected_map);
    const auto expected_values = expected.local_values_as_span();
    const auto actual_values = actual_on_expected_map.local_values_as_span();
    ASSERT_EQ(expected_values.size(), actual_values.size());
    for (std::size_t i = 0; i < expected_values.size(); ++i)
    {
      EXPECT_NEAR(expected_values[i], actual_values[i], tolerance) << "local vector entry " << i;
    }
  }

  TEST(ReducedLungDistributedTreeLinearSolverTests, SerialAirwaysMatchSparseSolverOnTwoRanks)
  {
    DistributedTreeFixture fixture;
    const double current_time = fixture.params.dynamics.time_increment;

    fixture.sync_state_from_x();
    auto residual = fixture.assemble_residual(current_time);
    fixture.assemble_jacobian(current_time);
    auto tree_linearization = fixture.assemble_tree_linearization(current_time);
    const auto tree_metadata = fixture.build_tree_metadata();

    Core::LinAlg::Vector<double> sparse_delta(*fixture.row_map, true);
    Core::LinAlg::Vector<double> tree_delta(*fixture.row_map, true);

    SparseNewtonLinearSolver sparse_solver(SparseNewtonLinearSolverContext{
        .comm = MPI_COMM_WORLD,
        .linear_solver_parameters = fixture.solver_params,
        .solver_params_callback = fixture.solver_params_callback,
        .correction_map = *fixture.row_map,
    });
    const NewtonLinearSystemMetadata linear_system_metadata{.current_time = current_time,
        .time_step_size_dt = fixture.params.dynamics.time_increment,
        .nonlinear_iteration = 0};
    sparse_solver.solve(
        *fixture.sysmat, residual, *fixture.x, linear_system_metadata, sparse_delta);

    DistributedTreeNewtonLinearSolver tree_solver(DistributedTreeNewtonLinearSolverContext{
        .tree_metadata = tree_metadata,
        .locally_relevant_dof_map = *fixture.locally_relevant_dof_map,
    });
    tree_solver.set_tree_linearization(tree_linearization);
    tree_solver.solve(*fixture.sysmat, residual, *fixture.x, linear_system_metadata, tree_delta);

    expect_vectors_near(sparse_delta, tree_delta, 1.0e-9);
  }
}  // namespace
