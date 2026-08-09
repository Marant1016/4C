// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_fem_discretization.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_rebalance.hpp"
#include "4C_reduced_lung_boundary_conditions.hpp"
#include "4C_reduced_lung_helpers.hpp"
#include "4C_reduced_lung_junctions.hpp"
#include "4C_reduced_lung_newton_linear_solver.hpp"
#include "4C_reduced_lung_newton_solver.hpp"
#include "4C_reduced_lung_terminal_unit.hpp"
#include "4C_utils_function_manager.hpp"
#include "4C_utils_function_of_time.hpp"

#include <mpi.h>
#include <Teuchos_ParameterList.hpp>

#include <algorithm>
#include <any>
#include <cmath>
#include <map>
#include <memory>
#include <numbers>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{
  using namespace FourC;
  using namespace FourC::ReducedLung;

  using ElementType = ReducedLungParameters::LungTree::ElementType;
  using BoundaryType = ReducedLungParameters::BoundaryConditions::Type;
  using RheologyType =
      ReducedLungParameters::LungTree::TerminalUnits::RheologicalModel::RheologicalModelType;
  using ElasticityType =
      ReducedLungParameters::LungTree::TerminalUnits::ElasticityModel::ElasticityModelType;

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

  ReducedLungParameters::Dynamics make_dynamics(double dt, int steps)
  {
    return ReducedLungParameters::Dynamics{
        .time_increment = dt,
        .number_of_steps = steps,
        .restart_every = 1,
        .results_every = 1,
        .linear_solver = 1,
        .max_nonlinear_iterations = 10,
        .nonlinear_residual_tolerance = 1e-8,
        .nonlinear_increment_tolerance = 1e-10,
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
        Core::IO::InputField<ReducedLungParameters::LungTree::Airways::FlowModel::ResistanceType>(
            ReducedLungParameters::LungTree::Airways::FlowModel::ResistanceType::Linear);
    params.lung_tree.airways.flow_model.include_inertia = Core::IO::InputField<bool>(false);
    params.lung_tree.airways.wall_model_type =
        Core::IO::InputField<ReducedLungParameters::LungTree::Airways::WallModelType>(
            ReducedLungParameters::LungTree::Airways::WallModelType::Rigid);
  }

  void set_dummy_terminal_unit_model(ReducedLungParameters& params)
  {
    params.lung_tree.terminal_units.rheological_model.rheological_model_type =
        Core::IO::InputField<RheologyType>(RheologyType::KelvinVoigt);
    params.lung_tree.terminal_units.rheological_model.kelvin_voigt.viscosity_kelvin_voigt_eta =
        Core::IO::InputField<double>(0.0);
    params.lung_tree.terminal_units.elasticity_model.elasticity_model_type =
        Core::IO::InputField<ElasticityType>(ElasticityType::Linear);
    params.lung_tree.terminal_units.elasticity_model.linear.elasticity_e =
        Core::IO::InputField<double>(1.0);
  }

  ReducedLungParameters make_terminal_unit_parameters(
      double dt, int steps, double radius, ElasticityType elasticity_type)
  {
    ReducedLungParameters params{};
    set_common_air_properties(params);
    params.dynamics = make_dynamics(dt, steps);

    params.lung_tree.topology.num_nodes = 2;
    params.lung_tree.topology.num_elements = 1;
    params.lung_tree.topology.node_coordinates =
        Core::IO::InputField<std::vector<double>>(std::unordered_map<int, std::vector<double>>{
            {1, {0.0, 0.0, 0.0}},
            {2, {radius, 0.0, 0.0}},
        });
    params.lung_tree.topology.element_nodes = Core::IO::InputField<std::vector<int>>(
        std::unordered_map<int, std::vector<int>>{{1, {1, 2}}});
    params.lung_tree.element_type = Core::IO::InputField<ElementType>(ElementType::TerminalUnit);
    params.lung_tree.generation = Core::IO::InputField<int>(-1);

    params.lung_tree.terminal_units.rheological_model.rheological_model_type =
        Core::IO::InputField<RheologyType>(RheologyType::KelvinVoigt);
    params.lung_tree.terminal_units.rheological_model.kelvin_voigt.viscosity_kelvin_voigt_eta =
        Core::IO::InputField<double>(1.0);
    params.lung_tree.terminal_units.elasticity_model.elasticity_model_type =
        Core::IO::InputField<ElasticityType>(elasticity_type);
    params.lung_tree.terminal_units.elasticity_model.linear.elasticity_e =
        Core::IO::InputField<double>(1.0);
    params.lung_tree.terminal_units.elasticity_model.ogden.ogden_parameter_kappa =
        Core::IO::InputField<double>(1.0);
    params.lung_tree.terminal_units.elasticity_model.ogden.ogden_parameter_beta =
        Core::IO::InputField<double>(-8.0);

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

  ReducedLungParameters make_serial_airway_parameters(double dt, int steps)
  {
    ReducedLungParameters params{};
    set_common_air_properties(params);
    params.dynamics = make_dynamics(dt, steps);

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
    set_linear_rigid_airway_model(params, {{1, 1.0}, {2, 0.9}, {3, 0.8}});
    set_dummy_terminal_unit_model(params);

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

  ReducedLungParameters make_bifurcation_airway_parameters(double dt, int steps)
  {
    ReducedLungParameters params{};
    set_common_air_properties(params);
    params.dynamics = make_dynamics(dt, steps);

    params.lung_tree.topology.num_nodes = 4;
    params.lung_tree.topology.num_elements = 3;
    params.lung_tree.topology.node_coordinates =
        Core::IO::InputField<std::vector<double>>(std::unordered_map<int, std::vector<double>>{
            {1, {0.0, 0.0, 0.0}},
            {2, {1.0, 0.0, 0.0}},
            {3, {2.0, 1.0, 0.0}},
            {4, {2.0, -1.0, 0.0}},
        });
    params.lung_tree.topology.element_nodes =
        Core::IO::InputField<std::vector<int>>(std::unordered_map<int, std::vector<int>>{
            {1, {1, 2}},
            {2, {2, 3}},
            {3, {2, 4}},
        });
    params.lung_tree.element_type = Core::IO::InputField<ElementType>(ElementType::Airway);
    params.lung_tree.generation =
        Core::IO::InputField<int>(std::unordered_map<int, int>{{1, 0}, {2, 1}, {3, 1}});
    set_linear_rigid_airway_model(params, {{1, 1.0}, {2, 0.85}, {3, 0.7}});
    set_dummy_terminal_unit_model(params);

    params.boundary_conditions.num_conditions = 3;
    params.boundary_conditions.bc_type = Core::IO::InputField<BoundaryType>(BoundaryType::Pressure);
    params.boundary_conditions.node_id =
        Core::IO::InputField<int>(std::unordered_map<int, int>{{1, 1}, {2, 3}, {3, 4}});
    params.boundary_conditions.value_source =
        ReducedLungParameters::BoundaryConditions::ValueSource::bc_function_id;
    params.boundary_conditions.function_id =
        Core::IO::InputField<int>(std::unordered_map<int, int>{{1, 1}, {2, 2}, {3, 2}});

    return params;
  }

  struct ReducedLungSolverFixture
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
    Teuchos::ParameterList solver_params;
    std::function<const Teuchos::ParameterList&(int)> solver_params_callback;

    ReducedLungSolverFixture(const std::string& name, const ReducedLungParameters& input_params,
        const std::vector<std::string>& function_definitions)
        : params(input_params),
          function_manager(make_function_manager(function_definitions)),
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
      sysmat = std::make_unique<Core::LinAlg::SparseMatrix>(*row_map, *locally_relevant_dof_map, 3);

      solver_params.set("SOLVER", Core::LinearSolver::SolverType::UMFPACK);
      solver_params.set("NAME", "Reduced_lung_solver");
      solver_params_callback = [this](int) -> const Teuchos::ParameterList&
      { return solver_params; };

      assembly_pipeline = create_default_reduced_lung_assembly_pipeline(
          airways, terminal_units, connections, bifurcations, boundary_conditions);
    }

    std::unique_ptr<NoxSolver> create_nox_solver()
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
      };
      return std::make_unique<NoxSolver>(context);
    }

    std::unique_ptr<NewtonSolver> create_newton_solver()
    {
      auto linear_solver =
          std::make_shared<SparseNewtonLinearSolver>(SparseNewtonLinearSolverContext{
              .comm = MPI_COMM_WORLD,
              .linear_solver_parameters = solver_params,
              .solver_params_callback = solver_params_callback,
              .correction_map = *row_map,
          });
      const NewtonSolverContext context{
          .dynamics = params.dynamics,
          .linear_solver = linear_solver,
          .assembly_pipeline = assembly_pipeline,
          .dofs = *dofs,
          .locally_relevant_dofs = *locally_relevant_dofs,
          .x = *x,
          .jacobian = *sysmat,
      };
      return std::make_unique<NewtonSolver>(context);
    }

    double residual_norm(double current_time)
    {
      Core::LinAlg::Vector<double> residual(*row_map, true);
      residual.put_scalar(0.0);
      for (const auto& assemble_residual : assembly_pipeline.residual_assemblers)
      {
        assemble_residual(
            residual, *locally_relevant_dofs, current_time, params.dynamics.time_increment);
      }

      double norm = 0.0;
      residual.norm_2(&norm);
      return norm;
    }

    void advance_end_of_timestep()
    {
      TerminalUnits::end_of_timestep_routine(
          terminal_units, *locally_relevant_dofs, params.dynamics.time_increment);
      Airways::end_of_timestep_routine(
          airways, *locally_relevant_dofs, params.dynamics.time_increment);
    }
  };

  void expect_vectors_near(const Core::LinAlg::Vector<double>& expected,
      const Core::LinAlg::Vector<double>& actual, double tolerance)
  {
    ASSERT_EQ(expected.global_length(), actual.global_length());
    Core::LinAlg::Vector<double> actual_on_expected_map(expected.get_map(), true);
    Core::LinAlg::export_to(actual, actual_on_expected_map);
    ASSERT_EQ(expected.local_length(), actual_on_expected_map.local_length());

    const auto expected_values = expected.local_values_as_span();
    const auto actual_values = actual_on_expected_map.local_values_as_span();
    for (std::size_t i = 0; i < expected_values.size(); ++i)
    {
      EXPECT_NEAR(expected_values[i], actual_values[i], tolerance) << "local vector entry " << i;
    }
  }

  void expect_terminal_unit_volumes_near(const ReducedLungSolverFixture& expected,
      const ReducedLungSolverFixture& actual, double tolerance)
  {
    ASSERT_EQ(expected.terminal_units.models.size(), actual.terminal_units.models.size());
    for (std::size_t model_index = 0; model_index < expected.terminal_units.models.size();
        ++model_index)
    {
      const auto& expected_data = expected.terminal_units.models[model_index].data;
      const auto& actual_data = actual.terminal_units.models[model_index].data;
      ASSERT_EQ(expected_data.volume_v.size(), actual_data.volume_v.size());
      for (std::size_t i = 0; i < expected_data.volume_v.size(); ++i)
      {
        EXPECT_NEAR(expected_data.volume_v[i], actual_data.volume_v[i], tolerance)
            << "terminal-unit model " << model_index << ", entry " << i;
      }
    }
  }

  void expect_connection_flow_balance(const ReducedLungSolverFixture& fixture, double tolerance)
  {
    const auto dofs = fixture.locally_relevant_dofs->local_values_as_span();
    for (std::size_t i = 0; i < fixture.connections.size(); ++i)
    {
      const auto& local_dof_ids = fixture.connections.local_dof_ids[i];
      const double flow_balance = dofs[local_dof_ids[Junctions::ConnectionData::q_out_parent]] -
                                  dofs[local_dof_ids[Junctions::ConnectionData::q_in_child]];
      EXPECT_NEAR(flow_balance, 0.0, tolerance) << "connection " << i;
    }
  }

  void expect_bifurcation_flow_balance(const ReducedLungSolverFixture& fixture, double tolerance)
  {
    const auto dofs = fixture.locally_relevant_dofs->local_values_as_span();
    for (std::size_t i = 0; i < fixture.bifurcations.size(); ++i)
    {
      const auto& local_dof_ids = fixture.bifurcations.local_dof_ids[i];
      const double flow_balance = dofs[local_dof_ids[Junctions::BifurcationData::q_out_parent]] -
                                  dofs[local_dof_ids[Junctions::BifurcationData::q_in_child_1]] -
                                  dofs[local_dof_ids[Junctions::BifurcationData::q_in_child_2]];
      EXPECT_NEAR(flow_balance, 0.0, tolerance) << "bifurcation " << i;
    }
  }

  struct ComparisonChecks
  {
    bool terminal_unit_volumes = false;
    bool connection_flow_balance = false;
    bool bifurcation_flow_balance = false;
  };

  void run_solver_comparison(const std::string& case_name, const ReducedLungParameters& params,
      const std::vector<std::string>& function_definitions, ComparisonChecks checks)
  {
    ReducedLungSolverFixture nox_fixture(case_name + "_nox", params, function_definitions);
    ReducedLungSolverFixture newton_fixture(case_name + "_newton", params, function_definitions);

    auto nox_solver = nox_fixture.create_nox_solver();
    auto newton_solver = newton_fixture.create_newton_solver();

    const double solution_tolerance = 1e-6;
    const double residual_tolerance =
        std::max(1e-7, 100.0 * params.dynamics.nonlinear_residual_tolerance);
    int completed_steps = 0;

    for (int step = 1; step <= params.dynamics.number_of_steps; ++step)
    {
      const double current_time = step * params.dynamics.time_increment;
      const unsigned int nox_iterations = nox_solver->solve(current_time);
      const unsigned int newton_iterations = newton_solver->solve(current_time);

      EXPECT_LE(
          nox_iterations, static_cast<unsigned int>(params.dynamics.max_nonlinear_iterations));
      EXPECT_LE(
          newton_iterations, static_cast<unsigned int>(params.dynamics.max_nonlinear_iterations));

      expect_vectors_near(*nox_fixture.x, *newton_fixture.x, solution_tolerance);
      expect_vectors_near(*nox_fixture.dofs, *newton_fixture.dofs, solution_tolerance);
      EXPECT_LE(nox_fixture.residual_norm(current_time), residual_tolerance);
      EXPECT_LE(newton_fixture.residual_norm(current_time), residual_tolerance);

      if (checks.connection_flow_balance)
      {
        expect_connection_flow_balance(nox_fixture, solution_tolerance);
        expect_connection_flow_balance(newton_fixture, solution_tolerance);
      }
      if (checks.bifurcation_flow_balance)
      {
        expect_bifurcation_flow_balance(nox_fixture, solution_tolerance);
        expect_bifurcation_flow_balance(newton_fixture, solution_tolerance);
      }

      nox_fixture.advance_end_of_timestep();
      newton_fixture.advance_end_of_timestep();
      if (checks.terminal_unit_volumes)
      {
        expect_terminal_unit_volumes_near(nox_fixture, newton_fixture, solution_tolerance);
      }

      ++completed_steps;
    }

    EXPECT_EQ(completed_steps, params.dynamics.number_of_steps);
  }

  TEST(ReducedLungNewtonVsNoxTest, SingleTerminalUnitOgdenMatchesNox)
  {
    const double radius = std::cbrt(3.0 / (4.0 * std::numbers::pi));
    const auto params = make_terminal_unit_parameters(0.4, 5, radius, ElasticityType::Ogden);
    run_solver_comparison("single_tu_ogden", params, {"-1/8*((1/(1+t)) - (1/(1+t))^(-7)) + 1", "0"},
        {.terminal_unit_volumes = true});
  }

  TEST(ReducedLungNewtonVsNoxTest, SingleTerminalUnitLinearElasticityMatchesNox)
  {
    const double radius = std::cbrt(3.0 / (4.0 * std::numbers::pi));
    const auto params = make_terminal_unit_parameters(0.25, 4, radius, ElasticityType::Linear);
    run_solver_comparison(
        "single_tu_linear", params, {"0.5*t", "0"}, {.terminal_unit_volumes = true});
  }

  TEST(ReducedLungNewtonVsNoxTest, SerialRigidAirwaysMatchNox)
  {
    const auto params = make_serial_airway_parameters(0.5, 3);
    run_solver_comparison("serial_airways", params, {"t", "0"}, {.connection_flow_balance = true});
  }

  TEST(ReducedLungNewtonVsNoxTest, BifurcationRigidAirwaysMatchNox)
  {
    const auto params = make_bifurcation_airway_parameters(0.5, 3);
    run_solver_comparison(
        "bifurcation_airways", params, {"t", "0"}, {.bifurcation_flow_balance = true});
  }
}  // namespace
