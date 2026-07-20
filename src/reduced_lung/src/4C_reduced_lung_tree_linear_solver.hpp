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

FOUR_C_NAMESPACE_OPEN

namespace ReducedLung
{
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
    const ReducedLungTreeMetadata& tree_metadata_;
    double pivot_tolerance_;
    TreeNewtonLinearSolverCoefficientSource coefficient_source_;
    const TreeLinearization* tree_linearization_ = nullptr;
  };
}  // namespace ReducedLung

FOUR_C_NAMESPACE_CLOSE

#endif
