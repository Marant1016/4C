// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_REDUCED_LUNG_TREE_LINEAR_SOLVER_HPP
#define FOUR_C_REDUCED_LUNG_TREE_LINEAR_SOLVER_HPP

#include "4C_config.hpp"

#include "4C_reduced_lung_linear_solver.hpp"
#include "4C_reduced_lung_tree_metadata.hpp"

#include <array>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace ReducedLung
{
  struct TreeNewtonLinearSolverProfile;

  struct TreeCoefficientLocation
  {
    int local_row = -1;
    int local_dof = -1;
    int structured_entry_index = -1;
  };

  enum class TreeNewtonLinearSolverCoefficientSource
  {
    SparseJacobian,
    StructuredTreeBlocks,
  };

  /**
   * @brief Context for the serial tree-based Newton correction solver.
   */
  struct TreeNewtonLinearSolverContext
  {
    const ReducedLungTreeMetadata& tree_metadata;
    double pivot_tolerance = 1.0e-12;
    TreeNewtonLinearSolverCoefficientSource coefficient_source =
        TreeNewtonLinearSolverCoefficientSource::SparseJacobian;
    TreeNewtonLinearSolverProfile* profile = nullptr;
  };

  /**
   * @brief Context for the distributed structured-tree Newton correction solver.
   */
  struct DistributedTreeNewtonLinearSolverContext
  {
    const ReducedLungTreeMetadata& tree_metadata;
    const Core::LinAlg::Map& locally_relevant_dof_map;
    double pivot_tolerance = 1.0e-12;
    TreeNewtonLinearSolverProfile* profile = nullptr;
  };

  /**
   * @brief Serial tree-based solver for reduced-lung Newton correction systems.
   *
   * The first implementation consumes the already assembled sparse Jacobian and residual, then uses
   * tree metadata to condense subtrees bottom-up and recover the correction top-down.
   */
  class TreeNewtonLinearSolver : public NewtonLinearSolver
  {
   public:
    explicit TreeNewtonLinearSolver(const TreeNewtonLinearSolverContext& context);

    [[nodiscard]] NewtonLinearizationType linearization_type() const override;

    void set_tree_linearization(const TreeLinearization& tree_linearization) override;

    void solve(Core::LinAlg::SparseMatrix& jacobian, const Core::LinAlg::Vector<double>& residual,
        const Core::LinAlg::Vector<double>& x, const NewtonLinearSystemMetadata& metadata,
        Core::LinAlg::Vector<double>& delta) override;

   private:
    struct ElementGroup
    {
      int begin = 0;
      int end = 0;
      int block_size = 0;
      int child_count = 0;
    };

    void build_symbolic_plan();

    void resolve_structured_coefficient_locations(const TreeLinearization& tree_linearization);

    const ReducedLungTreeMetadata& tree_metadata_;
    double pivot_tolerance_;
    TreeNewtonLinearSolverCoefficientSource coefficient_source_;
    const TreeLinearization* tree_linearization_ = nullptr;
    TreeNewtonLinearSolverProfile* profile_ = nullptr;

    int root_boundary_row_ = -1;
    int root_inlet_pressure_local_dof_ = -1;
    std::vector<int> global_element_id_;
    std::vector<int> inlet_pressure_local_dof_;
    std::vector<int> inlet_flow_unknown_index_;
    std::vector<int> outlet_pressure_unknown_index_;
    std::vector<int> block_size_;
    std::vector<int> child_interface_count_;
    std::vector<unsigned char> is_leaf_;

    std::vector<int> unknown_offset_;
    std::vector<int> equation_offset_;
    std::vector<int> child_interface_offset_;
    std::vector<int> matrix_offset_;

    std::vector<int> unknown_global_dof_ids_;
    std::vector<int> unknown_local_dof_ids_;
    std::vector<int> equation_rows_;
    std::vector<int> inlet_pressure_correction_local_dof_ids_;
    std::vector<int> unknown_correction_local_dof_ids_;
    bool correction_local_dof_ids_initialized_ = false;

    std::vector<int> child_element_index_;
    std::vector<int> pressure_row_;
    std::vector<int> parent_outlet_pressure_local_dof_;
    std::vector<int> child_inlet_pressure_local_dof_;
    std::vector<int> child_inlet_flow_local_dof_;
    std::vector<int> parent_outlet_pressure_unknown_index_;

    TreeCoefficientLocation root_boundary_coefficient_;
    std::vector<TreeCoefficientLocation> equation_inlet_pressure_coefficients_;
    std::vector<TreeCoefficientLocation> matrix_coefficients_;
    std::vector<TreeCoefficientLocation> child_pressure_parent_coefficients_;
    std::vector<TreeCoefficientLocation> child_pressure_child_coefficients_;
    std::vector<TreeCoefficientLocation> child_flow_coefficients_;

    double root_boundary_coefficient_value_ = 0.0;
    std::vector<double> equation_inlet_pressure_coefficient_values_;
    std::vector<double> matrix_coefficient_values_;
    std::vector<double> child_pressure_parent_coefficient_values_;
    std::vector<double> child_pressure_child_coefficient_values_;
    std::vector<double> child_flow_coefficient_values_;

    std::vector<int> grouped_element_indices_;
    std::vector<std::vector<ElementGroup>> bottom_up_layer_groups_;
    std::vector<std::vector<ElementGroup>> top_down_layer_groups_;
    bool use_scalar_tree_solve_ = false;

    std::vector<double> workspace_matrix_;
    std::vector<double> workspace_rhs_constant_;
    std::vector<double> workspace_rhs_inlet_pressure_;
    std::vector<double> workspace_intercept_;
    std::vector<double> workspace_slope_;
    std::vector<double> child_pressure_slope_;
    std::vector<double> child_pressure_intercept_;
    std::vector<double> batch_2x2_a00_;
    std::vector<double> batch_2x2_a01_;
    std::vector<double> batch_2x2_a10_;
    std::vector<double> batch_2x2_a11_;
    std::vector<double> batch_2x2_rhs_constant0_;
    std::vector<double> batch_2x2_rhs_constant1_;
    std::vector<double> batch_2x2_rhs_inlet_pressure0_;
    std::vector<double> batch_2x2_rhs_inlet_pressure1_;
    std::vector<double> batch_2x2_intercept0_;
    std::vector<double> batch_2x2_intercept1_;
    std::vector<double> batch_2x2_slope0_;
    std::vector<double> batch_2x2_slope1_;
    std::vector<int> batch_2x2_fallback_lanes_;
    std::vector<double> batch_3x3_a00_;
    std::vector<double> batch_3x3_a01_;
    std::vector<double> batch_3x3_a02_;
    std::vector<double> batch_3x3_a10_;
    std::vector<double> batch_3x3_a11_;
    std::vector<double> batch_3x3_a12_;
    std::vector<double> batch_3x3_a20_;
    std::vector<double> batch_3x3_a21_;
    std::vector<double> batch_3x3_a22_;
    std::vector<double> batch_3x3_rhs_constant0_;
    std::vector<double> batch_3x3_rhs_constant1_;
    std::vector<double> batch_3x3_rhs_constant2_;
    std::vector<double> batch_3x3_rhs_inlet_pressure0_;
    std::vector<double> batch_3x3_rhs_inlet_pressure1_;
    std::vector<double> batch_3x3_rhs_inlet_pressure2_;
    std::vector<double> batch_3x3_intercept0_;
    std::vector<double> batch_3x3_intercept1_;
    std::vector<double> batch_3x3_intercept2_;
    std::vector<double> batch_3x3_slope0_;
    std::vector<double> batch_3x3_slope1_;
    std::vector<double> batch_3x3_slope2_;
    std::vector<int> batch_3x3_fallback_lanes_;
    std::vector<double> top_down_inlet_pressure_;
    std::vector<double> top_down_unknown_values_;
    std::vector<double> top_down_outlet_pressure_;
    std::vector<double> top_down_child_pressure_;
    std::vector<std::string> element_context_;

    std::vector<double> subtree_relation_G_;
    std::vector<double> subtree_relation_h_;
    std::vector<double> inlet_pressure_by_element_;
    std::vector<int> inlet_pressure_stamp_;
    int current_solve_stamp_ = 0;
  };

  /**
   * @brief MPI-capable structured-tree Newton correction solver.
   *
   * This implementation performs bottom-up and top-down work on the ranks that own tree elements,
   * exchanging only condensed child-subtree relations, inlet-pressure corrections, and final
   * correction values across rank boundaries.
   */
  class DistributedTreeNewtonLinearSolver : public NewtonLinearSolver
  {
   public:
    explicit DistributedTreeNewtonLinearSolver(
        const DistributedTreeNewtonLinearSolverContext& context);

    [[nodiscard]] NewtonLinearizationType linearization_type() const override;

    void set_tree_linearization(const TreeLinearization& tree_linearization) override;

    void solve(Core::LinAlg::SparseMatrix& jacobian, const Core::LinAlg::Vector<double>& residual,
        const Core::LinAlg::Vector<double>& x, const NewtonLinearSystemMetadata& metadata,
        Core::LinAlg::Vector<double>& delta) override;

   private:
    struct SubtreeRelation
    {
      double G = 0.0;
      double h = 0.0;
    };

    struct ChildInterfacePlan
    {
      int child_element_index = -1;
      int pressure_global_row = -1;
      int parent_outlet_pressure_global_dof = -1;
      int child_inlet_pressure_global_dof = -1;
      int child_inlet_flow_global_dof = -1;
      int parent_outlet_pressure_unknown_index = -1;
    };

    struct ElementSolvePlan
    {
      int element_index = -1;
      int global_element_id = -1;
      int inlet_pressure_global_dof = -1;
      int inlet_flow_unknown_index = -1;
      bool is_leaf = false;
      std::vector<int> unknown_global_dof_ids;
      std::vector<int> equation_global_rows;
      std::vector<ChildInterfacePlan> child_interfaces;
      std::string context;
    };

    struct ElementWorkspace
    {
      std::vector<std::vector<double>> matrix;
      std::vector<double> rhs_constant;
      std::vector<double> rhs_inlet_pressure;
      std::vector<double> intercept;
      std::vector<double> slope;
      std::vector<double> child_pressure_slope;
      std::vector<double> child_pressure_intercept;
    };

    void build_symbolic_plan();

    MPI_Comm comm_ = MPI_COMM_NULL;
    const ReducedLungTreeMetadata& tree_metadata_;
    const Core::LinAlg::Map& locally_relevant_dof_map_;
    double pivot_tolerance_ = 1.0e-12;
    const TreeLinearization* tree_linearization_ = nullptr;
    TreeNewtonLinearSolverProfile* profile_ = nullptr;

    int root_boundary_global_row_ = -1;
    int root_inlet_pressure_global_dof_ = -1;
    std::vector<ElementSolvePlan> element_plans_;
    std::vector<ElementWorkspace> element_workspaces_;
    std::vector<SubtreeRelation> subtree_relations_;
    std::vector<double> inlet_pressure_by_element_;
  };
}  // namespace ReducedLung

FOUR_C_NAMESPACE_CLOSE

#endif
