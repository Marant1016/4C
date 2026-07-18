// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_reduced_lung_tree_linear_solver.hpp"

#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_exceptions.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace ReducedLung
{
  namespace
  {
    struct SubtreeRelation
    {
      double G = 0.0;
      double h = 0.0;
    };

    struct ElementRecoveryData
    {
      std::vector<int> unknown_global_dof_ids;
      std::vector<double> slope_by_inlet_pressure;
      std::vector<double> intercept;
    };

    struct ChildInterfaceEquation
    {
      int child_element_index = -1;
      int pressure_row = -1;
      int parent_outlet_pressure_local_dof = -1;
      int child_inlet_pressure_local_dof = -1;
      int child_inlet_flow_local_dof = -1;
    };

    const TreeJunctionMetadata* find_junction_for_parent(
        const ReducedLungTreeMetadata& tree_metadata, int parent_element_index)
    {
      for (const auto& junction : tree_metadata.junctions)
      {
        if (junction.parent_element_index == parent_element_index)
        {
          return &junction;
        }
      }
      return nullptr;
    }

    std::vector<const TreeBoundaryConditionMetadata*> outlet_boundaries_for_element(
        const ReducedLungTreeMetadata& tree_metadata, int element_index)
    {
      std::vector<const TreeBoundaryConditionMetadata*> boundaries;
      for (const auto& boundary : tree_metadata.boundary_conditions)
      {
        if (boundary.element_index == element_index && boundary.side == TreeBoundarySide::Outlet)
        {
          boundaries.push_back(&boundary);
        }
      }
      return boundaries;
    }

    const TreeBoundaryConditionMetadata& root_inlet_boundary(
        const ReducedLungTreeMetadata& tree_metadata)
    {
      const TreeBoundaryConditionMetadata* root_boundary = nullptr;
      for (const auto& boundary : tree_metadata.boundary_conditions)
      {
        if (boundary.element_index == tree_metadata.root_element_index &&
            boundary.side == TreeBoundarySide::Inlet)
        {
          FOUR_C_ASSERT_ALWAYS(root_boundary == nullptr,
              "TreeNewtonLinearSolver found multiple root inlet boundary conditions.");
          root_boundary = &boundary;
        }
      }
      FOUR_C_ASSERT_ALWAYS(
          root_boundary != nullptr, "TreeNewtonLinearSolver requires a root inlet boundary.");
      return *root_boundary;
    }

    double rhs_value(const Core::LinAlg::Vector<double>& residual, int local_row)
    {
      FOUR_C_ASSERT_ALWAYS(local_row >= 0 && local_row < residual.local_length(),
          "TreeNewtonLinearSolver row {} is not locally available.", local_row);
      return -residual.local_values_as_span()[static_cast<std::size_t>(local_row)];
    }

    double matrix_value(
        const Core::LinAlg::SparseMatrix& jacobian, int local_row, int local_col, double tolerance)
    {
      FOUR_C_ASSERT_ALWAYS(local_row >= 0 && local_row < jacobian.num_my_rows(),
          "TreeNewtonLinearSolver matrix row {} is not locally available.", local_row);
      FOUR_C_ASSERT_ALWAYS(local_col >= 0,
          "TreeNewtonLinearSolver matrix column is not locally available for row {}.", local_row);

      int n_entries = 0;
      double* values = nullptr;
      int* columns = nullptr;
      jacobian.extract_my_row_view(local_row, n_entries, values, columns);

      for (int i = 0; i < n_entries; ++i)
      {
        if (columns[i] == local_col)
        {
          return values[i];
        }
      }

      (void)tolerance;
      return 0.0;
    }

    double required_matrix_value(const Core::LinAlg::SparseMatrix& jacobian, int local_row,
        int local_col, double tolerance, const std::string& context)
    {
      const double value = matrix_value(jacobian, local_row, local_col, tolerance);
      FOUR_C_ASSERT_ALWAYS(std::abs(value) > tolerance,
          "TreeNewtonLinearSolver missing or near-zero matrix coefficient for {}.", context);
      return value;
    }

    std::vector<double> solve_dense_system(std::vector<std::vector<double>> matrix,
        std::vector<double> rhs, double pivot_tolerance, const std::string& context)
    {
      const int n = static_cast<int>(rhs.size());
      FOUR_C_ASSERT_ALWAYS(static_cast<int>(matrix.size()) == n,
          "TreeNewtonLinearSolver dense system has inconsistent size for {}.", context);

      for (int pivot_col = 0; pivot_col < n; ++pivot_col)
      {
        int pivot_row = pivot_col;
        double pivot_abs = std::abs(
            matrix[static_cast<std::size_t>(pivot_col)][static_cast<std::size_t>(pivot_col)]);
        for (int row = pivot_col + 1; row < n; ++row)
        {
          const double candidate_abs =
              std::abs(matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(pivot_col)]);
          if (candidate_abs > pivot_abs)
          {
            pivot_abs = candidate_abs;
            pivot_row = row;
          }
        }

        FOUR_C_ASSERT_ALWAYS(pivot_abs > pivot_tolerance,
            "TreeNewtonLinearSolver found a singular or underconstrained local block for {}.",
            context);

        if (pivot_row != pivot_col)
        {
          std::swap(matrix[static_cast<std::size_t>(pivot_col)],
              matrix[static_cast<std::size_t>(pivot_row)]);
          std::swap(
              rhs[static_cast<std::size_t>(pivot_col)], rhs[static_cast<std::size_t>(pivot_row)]);
        }

        const double pivot =
            matrix[static_cast<std::size_t>(pivot_col)][static_cast<std::size_t>(pivot_col)];
        for (int row = pivot_col + 1; row < n; ++row)
        {
          const double factor =
              matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(pivot_col)] / pivot;
          if (std::abs(factor) <= std::numeric_limits<double>::epsilon())
          {
            continue;
          }
          matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(pivot_col)] = 0.0;
          for (int col = pivot_col + 1; col < n; ++col)
          {
            matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)] -=
                factor * matrix[static_cast<std::size_t>(pivot_col)][static_cast<std::size_t>(col)];
          }
          rhs[static_cast<std::size_t>(row)] -= factor * rhs[static_cast<std::size_t>(pivot_col)];
        }
      }

      std::vector<double> solution(static_cast<std::size_t>(n), 0.0);
      for (int row = n - 1; row >= 0; --row)
      {
        double value = rhs[static_cast<std::size_t>(row)];
        for (int col = row + 1; col < n; ++col)
        {
          value -= matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)] *
                   solution[static_cast<std::size_t>(col)];
        }
        solution[static_cast<std::size_t>(row)] =
            value / matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)];
      }

      return solution;
    }

    int unknown_index_for_global_dof(
        const std::vector<int>& unknown_global_dof_ids, int global_dof_id)
    {
      const auto it =
          std::find(unknown_global_dof_ids.begin(), unknown_global_dof_ids.end(), global_dof_id);
      FOUR_C_ASSERT_ALWAYS(it != unknown_global_dof_ids.end(),
          "TreeNewtonLinearSolver recovery data does not contain global dof {}.", global_dof_id);
      return static_cast<int>(std::distance(unknown_global_dof_ids.begin(), it));
    }

    ChildInterfaceEquation child_interface_equation(const ReducedLungTreeMetadata& tree_metadata,
        const TreeJunctionMetadata& junction, int child_slot)
    {
      const int child_element_index =
          junction.child_element_indices[static_cast<std::size_t>(child_slot)];
      const auto& parent =
          tree_metadata.elements[static_cast<std::size_t>(junction.parent_element_index)];
      const auto& child = tree_metadata.elements[static_cast<std::size_t>(child_element_index)];

      return ChildInterfaceEquation{
          .child_element_index = child_element_index,
          .pressure_row = junction.first_local_equation_id + child_slot,
          .parent_outlet_pressure_local_dof = parent.local_dof_ids[1],
          .child_inlet_pressure_local_dof = child.local_dof_ids[0],
          .child_inlet_flow_local_dof = child.local_dof_ids[2],
      };
    }

    void add_equation_row(std::vector<std::vector<double>>& local_matrix,
        std::vector<double>& constant_rhs, std::vector<double>& inlet_pressure_rhs,
        const std::vector<int>& unknown_global_dof_ids, const TreeElementMetadata& element,
        const Core::LinAlg::SparseMatrix& jacobian, const Core::LinAlg::Vector<double>& residual,
        int local_row, double rhs_shift, double pivot_tolerance)
    {
      const std::size_t row = local_matrix.size();
      local_matrix.emplace_back(unknown_global_dof_ids.size(), 0.0);
      constant_rhs.push_back(rhs_value(residual, local_row) - rhs_shift);
      inlet_pressure_rhs.push_back(
          -matrix_value(jacobian, local_row, element.local_dof_ids[0], pivot_tolerance));

      for (std::size_t i = 0; i < unknown_global_dof_ids.size(); ++i)
      {
        const int global_dof = unknown_global_dof_ids[i];
        const int local_dof =
            element.local_dof_ids[static_cast<std::size_t>(global_dof - element.first_global_dof)];
        local_matrix[row][i] = matrix_value(jacobian, local_row, local_dof, pivot_tolerance);
      }
    }

    void add_downstream_flow_equation(std::vector<std::vector<double>>& local_matrix,
        std::vector<double>& constant_rhs, std::vector<double>& inlet_pressure_rhs,
        const std::vector<int>& unknown_global_dof_ids, const TreeElementMetadata& element,
        const TreeJunctionMetadata& junction, const ReducedLungTreeMetadata& tree_metadata,
        const std::vector<SubtreeRelation>& subtree_relations,
        const Core::LinAlg::SparseMatrix& jacobian, const Core::LinAlg::Vector<double>& residual,
        double pivot_tolerance)
    {
      const int flow_row = junction.first_local_equation_id + junction.child_count;
      add_equation_row(local_matrix, constant_rhs, inlet_pressure_rhs, unknown_global_dof_ids,
          element, jacobian, residual, flow_row, 0.0, pivot_tolerance);

      auto& row = local_matrix.back();
      double rhs_shift = 0.0;
      for (int child_slot = 0; child_slot < junction.child_count; ++child_slot)
      {
        const auto child_interface = child_interface_equation(tree_metadata, junction, child_slot);
        const auto& child_relation =
            subtree_relations[static_cast<std::size_t>(child_interface.child_element_index)];

        const double pressure_parent_coeff = required_matrix_value(jacobian,
            child_interface.pressure_row, child_interface.parent_outlet_pressure_local_dof,
            pivot_tolerance, "pressure-continuity parent pressure");
        const double pressure_child_coeff = required_matrix_value(jacobian,
            child_interface.pressure_row, child_interface.child_inlet_pressure_local_dof,
            pivot_tolerance, "pressure-continuity child pressure");
        const double pressure_rhs = rhs_value(residual, child_interface.pressure_row);

        const double child_pressure_slope = -pressure_parent_coeff / pressure_child_coeff;
        const double child_pressure_intercept = pressure_rhs / pressure_child_coeff;
        const double child_flow_slope = child_relation.G * child_pressure_slope;
        const double child_flow_intercept =
            child_relation.G * child_pressure_intercept + child_relation.h;

        const double flow_child_coeff =
            required_matrix_value(jacobian, flow_row, child_interface.child_inlet_flow_local_dof,
                pivot_tolerance, "junction flow child-flow coefficient");

        const int parent_outlet_pressure_global_dof = element.global_dof_ids[1];
        const int parent_outlet_pressure_unknown_index =
            unknown_index_for_global_dof(unknown_global_dof_ids, parent_outlet_pressure_global_dof);
        row[static_cast<std::size_t>(parent_outlet_pressure_unknown_index)] +=
            flow_child_coeff * child_flow_slope;
        rhs_shift += flow_child_coeff * child_flow_intercept;
      }

      constant_rhs.back() -= rhs_shift;
    }

    void condense_element(const ReducedLungTreeMetadata& tree_metadata, int element_index,
        const Core::LinAlg::SparseMatrix& jacobian, const Core::LinAlg::Vector<double>& residual,
        double pivot_tolerance, std::vector<SubtreeRelation>& subtree_relations,
        std::vector<ElementRecoveryData>& recovery_data)
    {
      const auto& element = tree_metadata.elements[static_cast<std::size_t>(element_index)];

      std::vector<int> unknown_global_dof_ids;
      unknown_global_dof_ids.reserve(static_cast<std::size_t>(element.num_dofs - 1));
      for (int i = 1; i < element.num_dofs; ++i)
      {
        unknown_global_dof_ids.push_back(element.global_dof_ids[static_cast<std::size_t>(i)]);
      }

      std::vector<std::vector<double>> local_matrix;
      std::vector<double> constant_rhs;
      std::vector<double> inlet_pressure_rhs;

      for (int row_offset = 0; row_offset < element.num_state_equations; ++row_offset)
      {
        add_equation_row(local_matrix, constant_rhs, inlet_pressure_rhs, unknown_global_dof_ids,
            element, jacobian, residual, element.first_local_state_equation_id + row_offset, 0.0,
            pivot_tolerance);
      }

      if (element.is_leaf())
      {
        const auto outlet_boundaries = outlet_boundaries_for_element(tree_metadata, element_index);
        FOUR_C_ASSERT_ALWAYS(outlet_boundaries.size() == 1u,
            "TreeNewtonLinearSolver requires exactly one outlet boundary for leaf element {}.",
            element.global_element_id + 1);
        add_equation_row(local_matrix, constant_rhs, inlet_pressure_rhs, unknown_global_dof_ids,
            element, jacobian, residual, outlet_boundaries.front()->local_equation_id, 0.0,
            pivot_tolerance);
      }
      else
      {
        const TreeJunctionMetadata* junction =
            find_junction_for_parent(tree_metadata, element_index);
        FOUR_C_ASSERT_ALWAYS(junction != nullptr,
            "TreeNewtonLinearSolver found no junction metadata for parent element {}.",
            element.global_element_id + 1);
        add_downstream_flow_equation(local_matrix, constant_rhs, inlet_pressure_rhs,
            unknown_global_dof_ids, element, *junction, tree_metadata, subtree_relations, jacobian,
            residual, pivot_tolerance);
      }

      FOUR_C_ASSERT_ALWAYS(local_matrix.size() == unknown_global_dof_ids.size(),
          "TreeNewtonLinearSolver local block for element {} has {} equations for {} unknowns.",
          element.global_element_id + 1, local_matrix.size(), unknown_global_dof_ids.size());

      const std::string context = "element " + std::to_string(element.global_element_id + 1);
      const auto intercept =
          solve_dense_system(local_matrix, constant_rhs, pivot_tolerance, context);
      const auto slope =
          solve_dense_system(local_matrix, inlet_pressure_rhs, pivot_tolerance, context);

      const int inlet_flow_global_dof = element.global_dof_ids[2];
      const int inlet_flow_index =
          unknown_index_for_global_dof(unknown_global_dof_ids, inlet_flow_global_dof);
      subtree_relations[static_cast<std::size_t>(element_index)] = SubtreeRelation{
          .G = slope[static_cast<std::size_t>(inlet_flow_index)],
          .h = intercept[static_cast<std::size_t>(inlet_flow_index)],
      };
      recovery_data[static_cast<std::size_t>(element_index)] = ElementRecoveryData{
          .unknown_global_dof_ids = std::move(unknown_global_dof_ids),
          .slope_by_inlet_pressure = slope,
          .intercept = intercept,
      };
    }

    void set_delta_value(Core::LinAlg::Vector<double>& delta, int global_dof_id, double value)
    {
      delta.replace_global_value(global_dof_id, value);
    }

    double recovered_dof_value(const TreeElementMetadata& element,
        const ElementRecoveryData& recovery_data, int global_dof_id, double inlet_pressure)
    {
      if (global_dof_id == element.global_dof_ids[0])
      {
        return inlet_pressure;
      }
      const int unknown_index =
          unknown_index_for_global_dof(recovery_data.unknown_global_dof_ids, global_dof_id);
      return recovery_data.slope_by_inlet_pressure[static_cast<std::size_t>(unknown_index)] *
                 inlet_pressure +
             recovery_data.intercept[static_cast<std::size_t>(unknown_index)];
    }

    void recover_element_and_children(const ReducedLungTreeMetadata& tree_metadata,
        int element_index, double inlet_pressure, const ElementRecoveryData& recovery_data,
        const Core::LinAlg::SparseMatrix& jacobian, const Core::LinAlg::Vector<double>& residual,
        double pivot_tolerance, Core::LinAlg::Vector<double>& delta,
        std::vector<double>& inlet_pressure_by_element)
    {
      const auto& element = tree_metadata.elements[static_cast<std::size_t>(element_index)];
      set_delta_value(delta, element.global_dof_ids[0], inlet_pressure);
      for (std::size_t i = 0; i < recovery_data.unknown_global_dof_ids.size(); ++i)
      {
        const double value =
            recovery_data.slope_by_inlet_pressure[i] * inlet_pressure + recovery_data.intercept[i];
        set_delta_value(delta, recovery_data.unknown_global_dof_ids[i], value);
      }

      if (element.is_leaf())
      {
        return;
      }

      const TreeJunctionMetadata* junction = find_junction_for_parent(tree_metadata, element_index);
      FOUR_C_ASSERT_ALWAYS(junction != nullptr,
          "TreeNewtonLinearSolver found no junction while recovering element {}.",
          element.global_element_id + 1);

      const double outlet_pressure =
          recovered_dof_value(element, recovery_data, element.global_dof_ids[1], inlet_pressure);
      for (int child_slot = 0; child_slot < junction->child_count; ++child_slot)
      {
        const auto child_interface = child_interface_equation(tree_metadata, *junction, child_slot);
        const double pressure_parent_coeff = required_matrix_value(jacobian,
            child_interface.pressure_row, child_interface.parent_outlet_pressure_local_dof,
            pivot_tolerance, "top-down pressure-continuity parent pressure");
        const double pressure_child_coeff = required_matrix_value(jacobian,
            child_interface.pressure_row, child_interface.child_inlet_pressure_local_dof,
            pivot_tolerance, "top-down pressure-continuity child pressure");
        const double pressure_rhs = rhs_value(residual, child_interface.pressure_row);

        inlet_pressure_by_element[static_cast<std::size_t>(child_interface.child_element_index)] =
            (pressure_rhs - pressure_parent_coeff * outlet_pressure) / pressure_child_coeff;
      }
    }

    double solve_root_inlet_pressure(const ReducedLungTreeMetadata& tree_metadata,
        const Core::LinAlg::SparseMatrix& jacobian, const Core::LinAlg::Vector<double>& residual,
        double pivot_tolerance)
    {
      const auto& root_element =
          tree_metadata.elements[static_cast<std::size_t>(tree_metadata.root_element_index)];
      const auto& boundary = root_inlet_boundary(tree_metadata);
      const double root_boundary_coeff = required_matrix_value(jacobian, boundary.local_equation_id,
          root_element.local_dof_ids[0], pivot_tolerance, "root inlet boundary");
      return rhs_value(residual, boundary.local_equation_id) / root_boundary_coeff;
    }
  }  // namespace

  TreeNewtonLinearSolver::TreeNewtonLinearSolver(const TreeNewtonLinearSolverContext& context)
      : tree_metadata_(context.tree_metadata), pivot_tolerance_(context.pivot_tolerance)
  {
    FOUR_C_ASSERT_ALWAYS(pivot_tolerance_ > 0.0,
        "TreeNewtonLinearSolver requires a positive pivot tolerance, got {}.", pivot_tolerance_);
  }

  void TreeNewtonLinearSolver::solve(Core::LinAlg::SparseMatrix& jacobian,
      const Core::LinAlg::Vector<double>& residual, const Core::LinAlg::Vector<double>& x,
      const NewtonLinearSystemMetadata& metadata, Core::LinAlg::Vector<double>& delta)
  {
    (void)x;
    (void)metadata;

    int comm_size = 1;
    MPI_Comm_size(delta.get_comm(), &comm_size);
    FOUR_C_ASSERT_ALWAYS(comm_size == 1,
        "TreeNewtonLinearSolver currently supports only serial reduced-lung solves.");
    FOUR_C_ASSERT_ALWAYS(jacobian.filled(),
        "TreeNewtonLinearSolver requires a completed sparse Jacobian before solving.");
    FOUR_C_ASSERT_ALWAYS(residual.local_length() == tree_metadata_.num_global_equations,
        "TreeNewtonLinearSolver requires all residual rows to be locally available.");
    FOUR_C_ASSERT_ALWAYS(delta.local_length() == tree_metadata_.num_global_dofs,
        "TreeNewtonLinearSolver requires all correction dofs to be locally available.");

    delta.put_scalar(0.0);

    std::vector<SubtreeRelation> subtree_relations(tree_metadata_.elements.size());
    std::vector<ElementRecoveryData> recovery_data(tree_metadata_.elements.size());

    for (const auto& layer : tree_metadata_.bottom_up_layers)
    {
      for (const int element_index : layer)
      {
        condense_element(tree_metadata_, element_index, jacobian, residual, pivot_tolerance_,
            subtree_relations, recovery_data);
      }
    }

    std::vector<double> inlet_pressure_by_element(
        tree_metadata_.elements.size(), std::numeric_limits<double>::quiet_NaN());
    inlet_pressure_by_element[static_cast<std::size_t>(tree_metadata_.root_element_index)] =
        solve_root_inlet_pressure(tree_metadata_, jacobian, residual, pivot_tolerance_);

    for (const auto& layer : tree_metadata_.top_down_layers)
    {
      for (const int element_index : layer)
      {
        const double inlet_pressure =
            inlet_pressure_by_element[static_cast<std::size_t>(element_index)];
        FOUR_C_ASSERT_ALWAYS(!std::isnan(inlet_pressure),
            "TreeNewtonLinearSolver missing inlet-pressure correction for element {}.",
            tree_metadata_.elements[static_cast<std::size_t>(element_index)].global_element_id + 1);
        recover_element_and_children(tree_metadata_, element_index, inlet_pressure,
            recovery_data[static_cast<std::size_t>(element_index)], jacobian, residual,
            pivot_tolerance_, delta, inlet_pressure_by_element);
      }
    }
  }
}  // namespace ReducedLung

FOUR_C_NAMESPACE_CLOSE
