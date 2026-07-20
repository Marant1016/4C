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
#include "4C_reduced_lung_tree_linearization.hpp"
#include "4C_utils_exceptions.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace ReducedLung
{
  namespace
  {
    class TreeCoefficientProvider
    {
     public:
      virtual ~TreeCoefficientProvider() = default;

      [[nodiscard]] virtual double value(int local_row, int local_col, double tolerance) const = 0;
    };

    class SparseTreeCoefficientProvider : public TreeCoefficientProvider
    {
     public:
      explicit SparseTreeCoefficientProvider(const Core::LinAlg::SparseMatrix& jacobian)
          : jacobian_(jacobian)
      {
      }

      [[nodiscard]] double value(int local_row, int local_col, double tolerance) const override
      {
        FOUR_C_ASSERT_ALWAYS(local_row >= 0 && local_row < jacobian_.num_my_rows(),
            "TreeNewtonLinearSolver matrix row {} is not locally available.", local_row);
        FOUR_C_ASSERT_ALWAYS(local_col >= 0,
            "TreeNewtonLinearSolver matrix column is not locally available for row {}.", local_row);

        int n_entries = 0;
        double* values = nullptr;
        int* columns = nullptr;
        jacobian_.extract_my_row_view(local_row, n_entries, values, columns);

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

     private:
      const Core::LinAlg::SparseMatrix& jacobian_;
    };

    class StructuredTreeCoefficientProvider : public TreeCoefficientProvider
    {
     public:
      explicit StructuredTreeCoefficientProvider(const TreeLinearization& linearization)
          : linearization_(linearization)
      {
      }

      [[nodiscard]] double value(int local_row, int local_col, double tolerance) const override
      {
        (void)tolerance;
        return linearization_.value(local_row, local_col);
      }

     private:
      const TreeLinearization& linearization_;
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
        const TreeCoefficientProvider& coefficients, int local_row, int local_col, double tolerance)
    {
      return coefficients.value(local_row, local_col, tolerance);
    }

    double required_matrix_value(const TreeCoefficientProvider& coefficients, int local_row,
        int local_col, double tolerance, const std::string& context)
    {
      const double value = matrix_value(coefficients, local_row, local_col, tolerance);
      FOUR_C_ASSERT_ALWAYS(std::abs(value) > tolerance,
          "TreeNewtonLinearSolver missing or near-zero matrix coefficient for {}.", context);
      return value;
    }

    void solve_dense_system(std::vector<std::vector<double>>& matrix, std::vector<double>& rhs_a,
        std::vector<double>& rhs_b, std::vector<double>& solution_a,
        std::vector<double>& solution_b, double pivot_tolerance, const std::string& context)
    {
      const int n = static_cast<int>(rhs_a.size());
      FOUR_C_ASSERT_ALWAYS(static_cast<int>(rhs_b.size()) == n,
          "TreeNewtonLinearSolver dense system has inconsistent RHS size for {}.", context);
      FOUR_C_ASSERT_ALWAYS(
          static_cast<int>(solution_a.size()) == n && static_cast<int>(solution_b.size()) == n,
          "TreeNewtonLinearSolver dense system has inconsistent solution size for {}.", context);
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
          std::swap(rhs_a[static_cast<std::size_t>(pivot_col)],
              rhs_a[static_cast<std::size_t>(pivot_row)]);
          std::swap(rhs_b[static_cast<std::size_t>(pivot_col)],
              rhs_b[static_cast<std::size_t>(pivot_row)]);
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
          rhs_a[static_cast<std::size_t>(row)] -=
              factor * rhs_a[static_cast<std::size_t>(pivot_col)];
          rhs_b[static_cast<std::size_t>(row)] -=
              factor * rhs_b[static_cast<std::size_t>(pivot_col)];
        }
      }

      for (int row = n - 1; row >= 0; --row)
      {
        double value_a = rhs_a[static_cast<std::size_t>(row)];
        double value_b = rhs_b[static_cast<std::size_t>(row)];
        for (int col = row + 1; col < n; ++col)
        {
          const double matrix_entry =
              matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)];
          value_a -= matrix_entry * solution_a[static_cast<std::size_t>(col)];
          value_b -= matrix_entry * solution_b[static_cast<std::size_t>(col)];
        }
        const double diagonal =
            matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)];
        solution_a[static_cast<std::size_t>(row)] = value_a / diagonal;
        solution_b[static_cast<std::size_t>(row)] = value_b / diagonal;
      }
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

    void set_delta_value(Core::LinAlg::Vector<double>& delta, int global_dof_id, double value)
    {
      delta.replace_global_value(global_dof_id, value);
    }
  }  // namespace

  TreeNewtonLinearSolver::TreeNewtonLinearSolver(const TreeNewtonLinearSolverContext& context)
      : tree_metadata_(context.tree_metadata),
        pivot_tolerance_(context.pivot_tolerance),
        coefficient_source_(context.coefficient_source)
  {
    FOUR_C_ASSERT_ALWAYS(pivot_tolerance_ > 0.0,
        "TreeNewtonLinearSolver requires a positive pivot tolerance, got {}.", pivot_tolerance_);
    build_symbolic_plan();
  }

  void TreeNewtonLinearSolver::build_symbolic_plan()
  {
    const auto& root_element =
        tree_metadata_.elements[static_cast<std::size_t>(tree_metadata_.root_element_index)];
    const auto& root_boundary = root_inlet_boundary(tree_metadata_);
    root_boundary_row_ = root_boundary.local_equation_id;
    root_inlet_pressure_local_dof_ = root_element.local_dof_ids[0];

    element_plans_.clear();
    element_plans_.resize(tree_metadata_.elements.size());
    element_workspaces_.clear();
    element_workspaces_.resize(tree_metadata_.elements.size());
    subtree_relations_.assign(tree_metadata_.elements.size(), SubtreeRelation{});
    inlet_pressure_by_element_.assign(
        tree_metadata_.elements.size(), std::numeric_limits<double>::quiet_NaN());

    for (std::size_t element_index = 0; element_index < tree_metadata_.elements.size();
        ++element_index)
    {
      const auto& element = tree_metadata_.elements[element_index];
      auto& plan = element_plans_[element_index];
      plan.element_index = static_cast<int>(element_index);
      plan.global_element_id = element.global_element_id;
      plan.inlet_pressure_local_dof = element.local_dof_ids[0];
      plan.is_leaf = element.is_leaf();
      plan.context = "element " + std::to_string(element.global_element_id + 1);

      plan.unknown_global_dof_ids.clear();
      plan.unknown_local_dof_ids.clear();
      plan.unknown_global_dof_ids.reserve(static_cast<std::size_t>(element.num_dofs - 1));
      plan.unknown_local_dof_ids.reserve(static_cast<std::size_t>(element.num_dofs - 1));
      for (int i = 1; i < element.num_dofs; ++i)
      {
        plan.unknown_global_dof_ids.push_back(element.global_dof_ids[static_cast<std::size_t>(i)]);
        plan.unknown_local_dof_ids.push_back(element.local_dof_ids[static_cast<std::size_t>(i)]);
      }

      plan.inlet_flow_unknown_index =
          unknown_index_for_global_dof(plan.unknown_global_dof_ids, element.global_dof_ids[2]);

      plan.equation_rows.clear();
      plan.equation_rows.reserve(plan.unknown_global_dof_ids.size());
      for (int row_offset = 0; row_offset < element.num_state_equations; ++row_offset)
      {
        plan.equation_rows.push_back(element.first_local_state_equation_id + row_offset);
      }

      plan.child_interfaces.clear();
      if (element.is_leaf())
      {
        const auto outlet_boundaries =
            outlet_boundaries_for_element(tree_metadata_, static_cast<int>(element_index));
        FOUR_C_ASSERT_ALWAYS(outlet_boundaries.size() == 1u,
            "TreeNewtonLinearSolver requires exactly one outlet boundary for leaf element {}.",
            element.global_element_id + 1);
        plan.equation_rows.push_back(outlet_boundaries.front()->local_equation_id);
      }
      else
      {
        const TreeJunctionMetadata* junction =
            find_junction_for_parent(tree_metadata_, static_cast<int>(element_index));
        FOUR_C_ASSERT_ALWAYS(junction != nullptr,
            "TreeNewtonLinearSolver found no junction metadata for parent element {}.",
            element.global_element_id + 1);

        plan.equation_rows.push_back(junction->first_local_equation_id + junction->child_count);
        plan.child_interfaces.reserve(static_cast<std::size_t>(junction->child_count));
        for (int child_slot = 0; child_slot < junction->child_count; ++child_slot)
        {
          const int child_element_index =
              junction->child_element_indices[static_cast<std::size_t>(child_slot)];
          const auto& child =
              tree_metadata_.elements[static_cast<std::size_t>(child_element_index)];
          const int parent_outlet_pressure_global_dof = element.global_dof_ids[1];
          plan.child_interfaces.push_back(ChildInterfacePlan{
              .child_element_index = child_element_index,
              .pressure_row = junction->first_local_equation_id + child_slot,
              .parent_outlet_pressure_local_dof = element.local_dof_ids[1],
              .child_inlet_pressure_local_dof = child.local_dof_ids[0],
              .child_inlet_flow_local_dof = child.local_dof_ids[2],
              .parent_outlet_pressure_unknown_index = unknown_index_for_global_dof(
                  plan.unknown_global_dof_ids, parent_outlet_pressure_global_dof),
          });
        }
      }

      FOUR_C_ASSERT_ALWAYS(plan.equation_rows.size() == plan.unknown_global_dof_ids.size(),
          "TreeNewtonLinearSolver local block for element {} has {} equations for {} unknowns.",
          element.global_element_id + 1, plan.equation_rows.size(),
          plan.unknown_global_dof_ids.size());

      auto& workspace = element_workspaces_[element_index];
      const std::size_t block_size = plan.unknown_global_dof_ids.size();
      workspace.matrix.assign(block_size, std::vector<double>(block_size, 0.0));
      workspace.rhs_constant.assign(block_size, 0.0);
      workspace.rhs_inlet_pressure.assign(block_size, 0.0);
      workspace.intercept.assign(block_size, 0.0);
      workspace.slope.assign(block_size, 0.0);
      workspace.child_pressure_slope.assign(plan.child_interfaces.size(), 0.0);
      workspace.child_pressure_intercept.assign(plan.child_interfaces.size(), 0.0);
    }
  }

  NewtonLinearizationType TreeNewtonLinearSolver::linearization_type() const
  {
    if (coefficient_source_ == TreeNewtonLinearSolverCoefficientSource::StructuredTreeBlocks)
    {
      return NewtonLinearizationType::StructuredTreeBlocks;
    }
    return NewtonLinearizationType::SparseJacobian;
  }

  void TreeNewtonLinearSolver::set_tree_linearization(const TreeLinearization& tree_linearization)
  {
    tree_linearization_ = &tree_linearization;
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
    if (coefficient_source_ == TreeNewtonLinearSolverCoefficientSource::SparseJacobian)
    {
      FOUR_C_ASSERT_ALWAYS(jacobian.filled(),
          "TreeNewtonLinearSolver requires a completed sparse Jacobian before solving.");
    }
    else
    {
      FOUR_C_ASSERT_ALWAYS(tree_linearization_ != nullptr,
          "TreeNewtonLinearSolver requires a structured tree linearization before solving.");
      FOUR_C_ASSERT_ALWAYS(tree_linearization_->num_rows() == tree_metadata_.num_global_equations,
          "TreeNewtonLinearSolver structured linearization row count does not match metadata.");
      FOUR_C_ASSERT_ALWAYS(
          tree_linearization_->num_dofs() == tree_metadata_.num_locally_relevant_dofs,
          "TreeNewtonLinearSolver structured linearization dof count does not match metadata.");
    }
    FOUR_C_ASSERT_ALWAYS(residual.local_length() == tree_metadata_.num_global_equations,
        "TreeNewtonLinearSolver requires all residual rows to be locally available.");
    FOUR_C_ASSERT_ALWAYS(delta.local_length() == tree_metadata_.num_global_dofs,
        "TreeNewtonLinearSolver requires all correction dofs to be locally available.");

    delta.put_scalar(0.0);

    SparseTreeCoefficientProvider sparse_coefficients(jacobian);
    std::optional<StructuredTreeCoefficientProvider> structured_coefficients;
    const TreeCoefficientProvider* coefficients = &sparse_coefficients;
    if (coefficient_source_ == TreeNewtonLinearSolverCoefficientSource::StructuredTreeBlocks)
    {
      structured_coefficients.emplace(*tree_linearization_);
      coefficients = &*structured_coefficients;
    }

    const auto add_equation_row = [&](const ElementSolvePlan& plan, ElementWorkspace& workspace,
                                      int equation_index, int local_row, double rhs_shift)
    {
      auto& matrix_row = workspace.matrix[static_cast<std::size_t>(equation_index)];
      workspace.rhs_constant[static_cast<std::size_t>(equation_index)] =
          rhs_value(residual, local_row) - rhs_shift;
      workspace.rhs_inlet_pressure[static_cast<std::size_t>(equation_index)] =
          -matrix_value(*coefficients, local_row, plan.inlet_pressure_local_dof, pivot_tolerance_);

      for (std::size_t i = 0; i < plan.unknown_local_dof_ids.size(); ++i)
      {
        matrix_row[i] =
            matrix_value(*coefficients, local_row, plan.unknown_local_dof_ids[i], pivot_tolerance_);
      }
    };

    std::fill(subtree_relations_.begin(), subtree_relations_.end(), SubtreeRelation{});

    for (const auto& layer : tree_metadata_.bottom_up_layers)
    {
      for (const int element_index : layer)
      {
        const auto& plan = element_plans_[static_cast<std::size_t>(element_index)];
        auto& workspace = element_workspaces_[static_cast<std::size_t>(element_index)];

        for (std::size_t equation_index = 0; equation_index < plan.equation_rows.size();
            ++equation_index)
        {
          const int local_row = plan.equation_rows[equation_index];
          add_equation_row(plan, workspace, static_cast<int>(equation_index), local_row, 0.0);

          const bool is_downstream_flow_equation =
              !plan.is_leaf && equation_index + 1 == plan.equation_rows.size();
          if (!is_downstream_flow_equation)
          {
            continue;
          }

          auto& matrix_row = workspace.matrix[equation_index];
          double rhs_shift = 0.0;
          for (std::size_t child_interface_index = 0;
              child_interface_index < plan.child_interfaces.size(); ++child_interface_index)
          {
            const auto& child_interface = plan.child_interfaces[child_interface_index];
            const auto& child_relation =
                subtree_relations_[static_cast<std::size_t>(child_interface.child_element_index)];

            const double pressure_parent_coeff = required_matrix_value(*coefficients,
                child_interface.pressure_row, child_interface.parent_outlet_pressure_local_dof,
                pivot_tolerance_, "pressure-continuity parent pressure");
            const double pressure_child_coeff = required_matrix_value(*coefficients,
                child_interface.pressure_row, child_interface.child_inlet_pressure_local_dof,
                pivot_tolerance_, "pressure-continuity child pressure");
            const double pressure_rhs = rhs_value(residual, child_interface.pressure_row);

            const double child_pressure_slope = -pressure_parent_coeff / pressure_child_coeff;
            const double child_pressure_intercept = pressure_rhs / pressure_child_coeff;
            workspace.child_pressure_slope[child_interface_index] = child_pressure_slope;
            workspace.child_pressure_intercept[child_interface_index] = child_pressure_intercept;

            const double child_flow_slope = child_relation.G * child_pressure_slope;
            const double child_flow_intercept =
                child_relation.G * child_pressure_intercept + child_relation.h;

            const double flow_child_coeff = required_matrix_value(*coefficients, local_row,
                child_interface.child_inlet_flow_local_dof, pivot_tolerance_,
                "junction flow child-flow coefficient");

            matrix_row[static_cast<std::size_t>(
                child_interface.parent_outlet_pressure_unknown_index)] +=
                flow_child_coeff * child_flow_slope;
            rhs_shift += flow_child_coeff * child_flow_intercept;
          }

          workspace.rhs_constant[equation_index] -= rhs_shift;
        }

        solve_dense_system(workspace.matrix, workspace.rhs_constant, workspace.rhs_inlet_pressure,
            workspace.intercept, workspace.slope, pivot_tolerance_, plan.context);

        subtree_relations_[static_cast<std::size_t>(element_index)] = SubtreeRelation{
            .G = workspace.slope[static_cast<std::size_t>(plan.inlet_flow_unknown_index)],
            .h = workspace.intercept[static_cast<std::size_t>(plan.inlet_flow_unknown_index)],
        };
      }
    }

    std::fill(inlet_pressure_by_element_.begin(), inlet_pressure_by_element_.end(),
        std::numeric_limits<double>::quiet_NaN());
    const double root_boundary_coeff = required_matrix_value(*coefficients, root_boundary_row_,
        root_inlet_pressure_local_dof_, pivot_tolerance_, "root inlet boundary");
    inlet_pressure_by_element_[static_cast<std::size_t>(tree_metadata_.root_element_index)] =
        rhs_value(residual, root_boundary_row_) / root_boundary_coeff;

    for (const auto& layer : tree_metadata_.top_down_layers)
    {
      for (const int element_index : layer)
      {
        const double inlet_pressure =
            inlet_pressure_by_element_[static_cast<std::size_t>(element_index)];
        FOUR_C_ASSERT_ALWAYS(!std::isnan(inlet_pressure),
            "TreeNewtonLinearSolver missing inlet-pressure correction for element {}.",
            tree_metadata_.elements[static_cast<std::size_t>(element_index)].global_element_id + 1);
        const auto& element = tree_metadata_.elements[static_cast<std::size_t>(element_index)];
        const auto& plan = element_plans_[static_cast<std::size_t>(element_index)];
        const auto& workspace = element_workspaces_[static_cast<std::size_t>(element_index)];

        set_delta_value(delta, element.global_dof_ids[0], inlet_pressure);
        for (std::size_t i = 0; i < plan.unknown_global_dof_ids.size(); ++i)
        {
          const double value = workspace.slope[i] * inlet_pressure + workspace.intercept[i];
          set_delta_value(delta, plan.unknown_global_dof_ids[i], value);
        }

        if (plan.is_leaf)
        {
          continue;
        }

        FOUR_C_ASSERT_ALWAYS(!plan.child_interfaces.empty(),
            "TreeNewtonLinearSolver found no children while recovering element {}.",
            plan.global_element_id + 1);
        const int outlet_pressure_index =
            plan.child_interfaces.front().parent_outlet_pressure_unknown_index;
        const double outlet_pressure =
            workspace.slope[static_cast<std::size_t>(outlet_pressure_index)] * inlet_pressure +
            workspace.intercept[static_cast<std::size_t>(outlet_pressure_index)];
        for (std::size_t child_interface_index = 0;
            child_interface_index < plan.child_interfaces.size(); ++child_interface_index)
        {
          const auto& child_interface = plan.child_interfaces[child_interface_index];
          inlet_pressure_by_element_[static_cast<std::size_t>(
              child_interface.child_element_index)] =
              workspace.child_pressure_slope[child_interface_index] * outlet_pressure +
              workspace.child_pressure_intercept[child_interface_index];
        }
      }
    }
  }
}  // namespace ReducedLung

FOUR_C_NAMESPACE_CLOSE
