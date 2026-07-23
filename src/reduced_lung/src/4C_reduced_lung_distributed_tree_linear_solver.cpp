// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_linalg_map.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_reduced_lung_solver_profile.hpp"
#include "4C_reduced_lung_tree_linear_solver.hpp"
#include "4C_reduced_lung_tree_linearization.hpp"
#include "4C_utils_exceptions.hpp"

#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace ReducedLung
{
  namespace
  {
    using Clock = std::chrono::steady_clock;

    double elapsed_seconds(const Clock::time_point start)
    {
      return std::chrono::duration<double>(Clock::now() - start).count();
    }

    std::vector<int> displacements_from_counts(const std::vector<int>& counts)
    {
      std::vector<int> displacements(counts.size(), 0);
      for (std::size_t i = 1; i < counts.size(); ++i)
      {
        displacements[i] = displacements[i - 1] + counts[i - 1];
      }
      return displacements;
    }

    std::vector<int> all_gather_counts(int local_count, MPI_Comm comm)
    {
      int comm_size = 1;
      MPI_Comm_size(comm, &comm_size);
      std::vector<int> counts(static_cast<std::size_t>(comm_size), 0);
      MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
      return counts;
    }

    std::vector<int> all_gatherv_int(
        const std::vector<int>& local_values, const std::vector<int>& counts, MPI_Comm comm)
    {
      const auto displacements = displacements_from_counts(counts);
      const int total_count = std::accumulate(counts.begin(), counts.end(), 0);
      std::vector<int> global_values(static_cast<std::size_t>(total_count), 0);
      MPI_Allgatherv(local_values.empty() ? nullptr : local_values.data(),
          static_cast<int>(local_values.size()), MPI_INT, global_values.data(), counts.data(),
          displacements.data(), MPI_INT, comm);
      return global_values;
    }

    std::vector<double> all_gatherv_double(
        const std::vector<double>& local_values, const std::vector<int>& counts, MPI_Comm comm)
    {
      const auto displacements = displacements_from_counts(counts);
      const int total_count = std::accumulate(counts.begin(), counts.end(), 0);
      std::vector<double> global_values(static_cast<std::size_t>(total_count), 0.0);
      MPI_Allgatherv(local_values.empty() ? nullptr : local_values.data(),
          static_cast<int>(local_values.size()), MPI_DOUBLE, global_values.data(), counts.data(),
          displacements.data(), MPI_DOUBLE, comm);
      return global_values;
    }

    std::vector<double> gather_global_residual(
        const Core::LinAlg::Vector<double>& residual, int num_global_equations, MPI_Comm comm)
    {
      std::vector<double> global_residual(static_cast<std::size_t>(num_global_equations), 0.0);
      const auto values = residual.local_values_as_span();
      const auto& row_map = residual.get_map();
      for (int local_row = 0; local_row < residual.local_length(); ++local_row)
      {
        const int global_row = row_map.gid(local_row);
        FOUR_C_ASSERT_ALWAYS(global_row >= 0 && global_row < num_global_equations,
            "DistributedTreeNewtonLinearSolver residual row {} is outside [0, {}).", global_row,
            num_global_equations);
        global_residual[static_cast<std::size_t>(global_row)] =
            values[static_cast<std::size_t>(local_row)];
      }
      MPI_Allreduce(
          MPI_IN_PLACE, global_residual.data(), num_global_equations, MPI_DOUBLE, MPI_SUM, comm);
      return global_residual;
    }

    struct GatheredCoefficients
    {
      std::vector<int> rows;
      std::vector<int> columns;
      std::vector<double> values;
    };

    GatheredCoefficients gather_global_coefficients(const TreeLinearization& linearization,
        const Core::LinAlg::Map& row_map, const Core::LinAlg::Map& locally_relevant_dof_map,
        MPI_Comm comm)
    {
      std::vector<int> local_rows;
      std::vector<int> local_columns;
      std::vector<double> local_values;

      for (int local_row = 0; local_row < linearization.num_rows(); ++local_row)
      {
        const int global_row = row_map.gid(local_row);
        FOUR_C_ASSERT_ALWAYS(global_row >= 0,
            "DistributedTreeNewtonLinearSolver found invalid local row {} in tree linearization.",
            local_row);
        for (const auto& [local_column, value] : linearization.entries(local_row))
        {
          const int global_column = locally_relevant_dof_map.gid(local_column);
          FOUR_C_ASSERT_ALWAYS(global_column >= 0,
              "DistributedTreeNewtonLinearSolver found invalid local dof {} in tree linearization.",
              local_column);
          local_rows.push_back(global_row);
          local_columns.push_back(global_column);
          local_values.push_back(value);
        }
      }

      const auto counts = all_gather_counts(static_cast<int>(local_values.size()), comm);
      return GatheredCoefficients{.rows = all_gatherv_int(local_rows, counts, comm),
          .columns = all_gatherv_int(local_columns, counts, comm),
          .values = all_gatherv_double(local_values, counts, comm)};
    }

    class GlobalTreeCoefficientProvider
    {
     public:
      GlobalTreeCoefficientProvider(int num_global_dofs, const GatheredCoefficients& coefficients,
          TreeNewtonLinearSolverProfile* profile)
          : num_global_dofs_(num_global_dofs), profile_(profile)
      {
        FOUR_C_ASSERT_ALWAYS(coefficients.rows.size() == coefficients.columns.size() &&
                                 coefficients.rows.size() == coefficients.values.size(),
            "DistributedTreeNewtonLinearSolver gathered inconsistent coefficient arrays.");
        values_.reserve(coefficients.values.size());
        for (std::size_t i = 0; i < coefficients.values.size(); ++i)
        {
          values_[key(coefficients.rows[i], coefficients.columns[i])] = coefficients.values[i];
        }
      }

      [[nodiscard]] double value(int global_row, int global_column, double tolerance) const
      {
        (void)tolerance;
        const auto lookup_start = profile_ != nullptr ? Clock::now() : Clock::time_point{};
        const auto entry = values_.find(key(global_row, global_column));
        if (profile_ != nullptr)
        {
          profile_->coefficient_lookup_time += elapsed_seconds(lookup_start);
          ++profile_->coefficient_lookup_count;
        }
        return entry == values_.end() ? 0.0 : entry->second;
      }

     private:
      [[nodiscard]] std::int64_t key(int global_row, int global_column) const
      {
        FOUR_C_ASSERT_ALWAYS(global_row >= 0 && global_column >= 0,
            "DistributedTreeNewtonLinearSolver coefficient ids must be non-negative.");
        return static_cast<std::int64_t>(global_row) * static_cast<std::int64_t>(num_global_dofs_) +
               static_cast<std::int64_t>(global_column);
      }

      int num_global_dofs_ = 0;
      TreeNewtonLinearSolverProfile* profile_ = nullptr;
      std::unordered_map<std::int64_t, double> values_;
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
              "DistributedTreeNewtonLinearSolver found multiple root inlet boundary conditions.");
          root_boundary = &boundary;
        }
      }
      FOUR_C_ASSERT_ALWAYS(root_boundary != nullptr,
          "DistributedTreeNewtonLinearSolver requires a root inlet boundary.");
      return *root_boundary;
    }

    double rhs_value(const std::vector<double>& global_residual, int global_row)
    {
      FOUR_C_ASSERT_ALWAYS(global_row >= 0 && global_row < static_cast<int>(global_residual.size()),
          "DistributedTreeNewtonLinearSolver row {} is not globally available.", global_row);
      return -global_residual[static_cast<std::size_t>(global_row)];
    }

    double required_matrix_value(const GlobalTreeCoefficientProvider& coefficients, int global_row,
        int global_column, double tolerance, const std::string& context)
    {
      const double value = coefficients.value(global_row, global_column, tolerance);
      FOUR_C_ASSERT_ALWAYS(std::abs(value) > tolerance,
          "DistributedTreeNewtonLinearSolver missing or near-zero matrix coefficient for {}.",
          context);
      return value;
    }

    void solve_dense_system(std::vector<std::vector<double>>& matrix, std::vector<double>& rhs_a,
        std::vector<double>& rhs_b, std::vector<double>& solution_a,
        std::vector<double>& solution_b, double pivot_tolerance, const std::string& context)
    {
      const int n = static_cast<int>(rhs_a.size());
      FOUR_C_ASSERT_ALWAYS(static_cast<int>(rhs_b.size()) == n,
          "DistributedTreeNewtonLinearSolver dense system has inconsistent RHS size for {}.",
          context);
      FOUR_C_ASSERT_ALWAYS(
          static_cast<int>(solution_a.size()) == n && static_cast<int>(solution_b.size()) == n,
          "DistributedTreeNewtonLinearSolver dense system has inconsistent solution size for {}.",
          context);
      FOUR_C_ASSERT_ALWAYS(static_cast<int>(matrix.size()) == n,
          "DistributedTreeNewtonLinearSolver dense system has inconsistent size for {}.", context);

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
            "DistributedTreeNewtonLinearSolver found a singular or underconstrained local block "
            "for {}.",
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
          "DistributedTreeNewtonLinearSolver recovery data does not contain global dof {}.",
          global_dof_id);
      return static_cast<int>(std::distance(unknown_global_dof_ids.begin(), it));
    }

    void set_local_delta_value(Core::LinAlg::Vector<double>& delta, int global_dof_id, double value)
    {
      if (delta.get_map().lid(global_dof_id) >= 0)
      {
        delta.replace_global_value(global_dof_id, value);
      }
    }
  }  // namespace

  DistributedTreeNewtonLinearSolver::DistributedTreeNewtonLinearSolver(
      const DistributedTreeNewtonLinearSolverContext& context)
      : comm_(context.locally_relevant_dof_map.get_comm()),
        tree_metadata_(context.tree_metadata),
        locally_relevant_dof_map_(context.locally_relevant_dof_map),
        pivot_tolerance_(context.pivot_tolerance),
        profile_(context.profile)
  {
    FOUR_C_ASSERT_ALWAYS(pivot_tolerance_ > 0.0,
        "DistributedTreeNewtonLinearSolver requires a positive pivot tolerance, got {}.",
        pivot_tolerance_);
    build_symbolic_plan();
  }

  void DistributedTreeNewtonLinearSolver::build_symbolic_plan()
  {
    const auto& root_element =
        tree_metadata_.elements[static_cast<std::size_t>(tree_metadata_.root_element_index)];
    const auto& root_boundary = root_inlet_boundary(tree_metadata_);
    root_boundary_global_row_ = root_boundary.global_equation_id;
    root_inlet_pressure_global_dof_ = root_element.global_dof_ids[0];

    element_plans_.clear();
    element_plans_.resize(tree_metadata_.elements.size());
    element_workspaces_.clear();
    element_workspaces_.resize(tree_metadata_.elements.size());
    subtree_relations_.assign(tree_metadata_.elements.size(), SubtreeRelation{});
    inlet_pressure_by_element_.assign(
        tree_metadata_.elements.size(), std::numeric_limits<double>::quiet_NaN());
    if (profile_ != nullptr)
    {
      profile_->element_count = 0;
      profile_->total_local_block_dofs = 0;
      profile_->max_local_block_size = 0;
    }

    for (std::size_t element_index = 0; element_index < tree_metadata_.elements.size();
        ++element_index)
    {
      const auto& element = tree_metadata_.elements[element_index];
      auto& plan = element_plans_[element_index];
      plan.element_index = static_cast<int>(element_index);
      plan.global_element_id = element.global_element_id;
      plan.inlet_pressure_global_dof = element.global_dof_ids[0];
      plan.is_leaf = element.is_leaf();
      plan.context = "element " + std::to_string(element.global_element_id + 1);

      plan.unknown_global_dof_ids.clear();
      plan.unknown_global_dof_ids.reserve(static_cast<std::size_t>(element.num_dofs - 1));
      for (int i = 1; i < element.num_dofs; ++i)
      {
        plan.unknown_global_dof_ids.push_back(element.global_dof_ids[static_cast<std::size_t>(i)]);
      }
      plan.inlet_flow_unknown_index =
          unknown_index_for_global_dof(plan.unknown_global_dof_ids, element.global_dof_ids[2]);

      plan.equation_global_rows.clear();
      plan.equation_global_rows.reserve(plan.unknown_global_dof_ids.size());
      for (int row_offset = 0; row_offset < element.num_state_equations; ++row_offset)
      {
        plan.equation_global_rows.push_back(element.first_global_state_equation_id + row_offset);
      }

      plan.child_interfaces.clear();
      if (element.is_leaf())
      {
        const auto outlet_boundaries =
            outlet_boundaries_for_element(tree_metadata_, static_cast<int>(element_index));
        FOUR_C_ASSERT_ALWAYS(outlet_boundaries.size() == 1u,
            "DistributedTreeNewtonLinearSolver requires exactly one outlet boundary for leaf "
            "element {}.",
            element.global_element_id + 1);
        plan.equation_global_rows.push_back(outlet_boundaries.front()->global_equation_id);
      }
      else
      {
        const TreeJunctionMetadata* junction =
            find_junction_for_parent(tree_metadata_, static_cast<int>(element_index));
        FOUR_C_ASSERT_ALWAYS(junction != nullptr,
            "DistributedTreeNewtonLinearSolver found no junction metadata for parent element {}.",
            element.global_element_id + 1);

        plan.equation_global_rows.push_back(
            junction->first_global_equation_id + junction->child_count);
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
              .pressure_global_row = junction->first_global_equation_id + child_slot,
              .parent_outlet_pressure_global_dof = parent_outlet_pressure_global_dof,
              .child_inlet_pressure_global_dof = child.global_dof_ids[0],
              .child_inlet_flow_global_dof = child.global_dof_ids[2],
              .parent_outlet_pressure_unknown_index = unknown_index_for_global_dof(
                  plan.unknown_global_dof_ids, parent_outlet_pressure_global_dof),
          });
        }
      }

      FOUR_C_ASSERT_ALWAYS(plan.equation_global_rows.size() == plan.unknown_global_dof_ids.size(),
          "DistributedTreeNewtonLinearSolver local block for element {} has {} equations for {} "
          "unknowns.",
          element.global_element_id + 1, plan.equation_global_rows.size(),
          plan.unknown_global_dof_ids.size());

      auto& workspace = element_workspaces_[element_index];
      const std::size_t block_size = plan.unknown_global_dof_ids.size();
      if (profile_ != nullptr)
      {
        ++profile_->element_count;
        profile_->total_local_block_dofs += block_size;
        profile_->max_local_block_size =
            std::max(profile_->max_local_block_size, static_cast<int>(block_size));
      }
      workspace.matrix.assign(block_size, std::vector<double>(block_size, 0.0));
      workspace.rhs_constant.assign(block_size, 0.0);
      workspace.rhs_inlet_pressure.assign(block_size, 0.0);
      workspace.intercept.assign(block_size, 0.0);
      workspace.slope.assign(block_size, 0.0);
      workspace.child_pressure_slope.assign(plan.child_interfaces.size(), 0.0);
      workspace.child_pressure_intercept.assign(plan.child_interfaces.size(), 0.0);
    }
  }

  NewtonLinearizationType DistributedTreeNewtonLinearSolver::linearization_type() const
  {
    return NewtonLinearizationType::StructuredTreeBlocks;
  }

  void DistributedTreeNewtonLinearSolver::set_tree_linearization(
      const TreeLinearization& tree_linearization)
  {
    tree_linearization_ = &tree_linearization;
  }

  void DistributedTreeNewtonLinearSolver::solve(Core::LinAlg::SparseMatrix& jacobian,
      const Core::LinAlg::Vector<double>& residual, const Core::LinAlg::Vector<double>& x,
      const NewtonLinearSystemMetadata& metadata, Core::LinAlg::Vector<double>& delta)
  {
    (void)jacobian;
    (void)x;
    (void)metadata;
    const auto solve_start = Clock::now();

    FOUR_C_ASSERT_ALWAYS(tree_linearization_ != nullptr,
        "DistributedTreeNewtonLinearSolver requires a structured tree linearization before "
        "solving.");
    FOUR_C_ASSERT_ALWAYS(tree_linearization_->num_rows() == residual.local_length(),
        "DistributedTreeNewtonLinearSolver structured linearization row count does not match the "
        "local residual.");
    FOUR_C_ASSERT_ALWAYS(
        tree_linearization_->num_dofs() == locally_relevant_dof_map_.num_my_elements(),
        "DistributedTreeNewtonLinearSolver structured linearization dof count does not match the "
        "local column map.");
    FOUR_C_ASSERT_ALWAYS(residual.global_length() == tree_metadata_.num_global_equations,
        "DistributedTreeNewtonLinearSolver residual size does not match tree metadata.");
    FOUR_C_ASSERT_ALWAYS(delta.global_length() == tree_metadata_.num_global_dofs,
        "DistributedTreeNewtonLinearSolver correction size does not match tree metadata.");

    delta.put_scalar(0.0);

    const auto communication_start = Clock::now();
    const auto global_residual =
        gather_global_residual(residual, tree_metadata_.num_global_equations, comm_);
    const auto gathered_coefficients = gather_global_coefficients(
        *tree_linearization_, residual.get_map(), locally_relevant_dof_map_, comm_);
    if (profile_ != nullptr)
    {
      profile_->communication_time += elapsed_seconds(communication_start);
      profile_->communicated_residual_count += static_cast<std::uint64_t>(residual.local_length());
      profile_->communicated_coefficient_count +=
          static_cast<std::uint64_t>(gathered_coefficients.values.size());
      profile_->communication_bytes += static_cast<std::uint64_t>(
          global_residual.size() * sizeof(double) +
          gathered_coefficients.values.size() * (2 * sizeof(int) + sizeof(double)));
    }
    const GlobalTreeCoefficientProvider coefficients(
        tree_metadata_.num_global_dofs, gathered_coefficients, profile_);

    const auto add_equation_row = [&](const ElementSolvePlan& plan, ElementWorkspace& workspace,
                                      int equation_index, int global_row, double rhs_shift)
    {
      auto& matrix_row = workspace.matrix[static_cast<std::size_t>(equation_index)];
      workspace.rhs_constant[static_cast<std::size_t>(equation_index)] =
          rhs_value(global_residual, global_row) - rhs_shift;
      workspace.rhs_inlet_pressure[static_cast<std::size_t>(equation_index)] =
          -coefficients.value(global_row, plan.inlet_pressure_global_dof, pivot_tolerance_);

      for (std::size_t i = 0; i < plan.unknown_global_dof_ids.size(); ++i)
      {
        matrix_row[i] =
            coefficients.value(global_row, plan.unknown_global_dof_ids[i], pivot_tolerance_);
      }
    };

    std::fill(subtree_relations_.begin(), subtree_relations_.end(), SubtreeRelation{});

    const auto bottom_up_start = Clock::now();
    for (const auto& layer : tree_metadata_.bottom_up_layers)
    {
      for (const int element_index : layer)
      {
        const auto& plan = element_plans_[static_cast<std::size_t>(element_index)];
        auto& workspace = element_workspaces_[static_cast<std::size_t>(element_index)];

        for (std::size_t equation_index = 0; equation_index < plan.equation_global_rows.size();
            ++equation_index)
        {
          const int global_row = plan.equation_global_rows[equation_index];
          add_equation_row(plan, workspace, static_cast<int>(equation_index), global_row, 0.0);

          const bool is_downstream_flow_equation =
              !plan.is_leaf && equation_index + 1 == plan.equation_global_rows.size();
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

            const double pressure_parent_coeff =
                required_matrix_value(coefficients, child_interface.pressure_global_row,
                    child_interface.parent_outlet_pressure_global_dof, pivot_tolerance_,
                    "pressure-continuity parent pressure");
            const double pressure_child_coeff =
                required_matrix_value(coefficients, child_interface.pressure_global_row,
                    child_interface.child_inlet_pressure_global_dof, pivot_tolerance_,
                    "pressure-continuity child pressure");
            const double pressure_rhs =
                rhs_value(global_residual, child_interface.pressure_global_row);

            const double child_pressure_slope = -pressure_parent_coeff / pressure_child_coeff;
            const double child_pressure_intercept = pressure_rhs / pressure_child_coeff;
            workspace.child_pressure_slope[child_interface_index] = child_pressure_slope;
            workspace.child_pressure_intercept[child_interface_index] = child_pressure_intercept;

            const double child_flow_slope = child_relation.G * child_pressure_slope;
            const double child_flow_intercept =
                child_relation.G * child_pressure_intercept + child_relation.h;

            const double flow_child_coeff = required_matrix_value(coefficients, global_row,
                child_interface.child_inlet_flow_global_dof, pivot_tolerance_,
                "junction flow child-flow coefficient");

            matrix_row[static_cast<std::size_t>(
                child_interface.parent_outlet_pressure_unknown_index)] +=
                flow_child_coeff * child_flow_slope;
            rhs_shift += flow_child_coeff * child_flow_intercept;
          }

          workspace.rhs_constant[equation_index] -= rhs_shift;
        }

        const auto dense_solve_start = Clock::now();
        solve_dense_system(workspace.matrix, workspace.rhs_constant, workspace.rhs_inlet_pressure,
            workspace.intercept, workspace.slope, pivot_tolerance_, plan.context);
        if (profile_ != nullptr)
        {
          profile_->dense_solve_time += elapsed_seconds(dense_solve_start);
          ++profile_->dense_solve_count;
        }

        subtree_relations_[static_cast<std::size_t>(element_index)] = SubtreeRelation{
            .G = workspace.slope[static_cast<std::size_t>(plan.inlet_flow_unknown_index)],
            .h = workspace.intercept[static_cast<std::size_t>(plan.inlet_flow_unknown_index)],
        };
      }
    }
    if (profile_ != nullptr)
    {
      profile_->bottom_up_time += elapsed_seconds(bottom_up_start);
    }

    std::fill(inlet_pressure_by_element_.begin(), inlet_pressure_by_element_.end(),
        std::numeric_limits<double>::quiet_NaN());
    const double root_boundary_coeff =
        required_matrix_value(coefficients, root_boundary_global_row_,
            root_inlet_pressure_global_dof_, pivot_tolerance_, "root inlet boundary");
    inlet_pressure_by_element_[static_cast<std::size_t>(tree_metadata_.root_element_index)] =
        rhs_value(global_residual, root_boundary_global_row_) / root_boundary_coeff;

    const auto top_down_start = Clock::now();
    for (const auto& layer : tree_metadata_.top_down_layers)
    {
      for (const int element_index : layer)
      {
        const double inlet_pressure =
            inlet_pressure_by_element_[static_cast<std::size_t>(element_index)];
        FOUR_C_ASSERT_ALWAYS(!std::isnan(inlet_pressure),
            "DistributedTreeNewtonLinearSolver missing inlet-pressure correction for element {}.",
            tree_metadata_.elements[static_cast<std::size_t>(element_index)].global_element_id + 1);
        const auto& element = tree_metadata_.elements[static_cast<std::size_t>(element_index)];
        const auto& plan = element_plans_[static_cast<std::size_t>(element_index)];
        const auto& workspace = element_workspaces_[static_cast<std::size_t>(element_index)];

        set_local_delta_value(delta, element.global_dof_ids[0], inlet_pressure);
        for (std::size_t i = 0; i < plan.unknown_global_dof_ids.size(); ++i)
        {
          const double value = workspace.slope[i] * inlet_pressure + workspace.intercept[i];
          set_local_delta_value(delta, plan.unknown_global_dof_ids[i], value);
        }

        if (plan.is_leaf)
        {
          continue;
        }

        FOUR_C_ASSERT_ALWAYS(!plan.child_interfaces.empty(),
            "DistributedTreeNewtonLinearSolver found no children while recovering element {}.",
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
    if (profile_ != nullptr)
    {
      profile_->top_down_time += elapsed_seconds(top_down_start);
      profile_->total_solve_time += elapsed_seconds(solve_start);
      ++profile_->solve_count;
    }
  }
}  // namespace ReducedLung

FOUR_C_NAMESPACE_CLOSE
