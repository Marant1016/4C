// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_REDUCED_LUNG_LINEAR_SOLVER_HPP
#define FOUR_C_REDUCED_LUNG_LINEAR_SOLVER_HPP

#include "4C_config.hpp"

#include <mpi.h>
#include <Teuchos_ParameterList.hpp>

#include <functional>
#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class Map;
  class Solver;
  class SparseMatrix;
  template <typename T>
  class Vector;
}  // namespace Core::LinAlg

namespace ReducedLung
{
  struct SparseNewtonLinearSolverProfile;
  class TreeCoefficientAssemblyTarget;
  class TreeLinearization;

  /**
   * @brief Metadata passed to linear solvers for reduced-lung Newton correction systems.
   */
  struct NewtonLinearSystemMetadata
  {
    double current_time = 0.0;
    double time_step_size_dt = 0.0;
    unsigned int nonlinear_iteration = 0;
  };

  enum class NewtonLinearizationType
  {
    SparseJacobian,
    StructuredTreeBlocks,
  };

  /**
   * @brief Interface for reduced-lung Newton correction linear solvers.
   */
  class NewtonLinearSolver
  {
   public:
    virtual ~NewtonLinearSolver() = default;

    [[nodiscard]] virtual NewtonLinearizationType linearization_type() const
    {
      return NewtonLinearizationType::SparseJacobian;
    }

    virtual void set_tree_linearization(const TreeLinearization& tree_linearization)
    {
      (void)tree_linearization;
    }

    [[nodiscard]] virtual TreeCoefficientAssemblyTarget* direct_tree_coefficient_target()
    {
      return nullptr;
    }

    /**
     * @brief Solve one Newton correction system.
     *
     * Sign convention: solve `jacobian * delta = -residual`. The nonlinear solver then applies
     * `x = x + delta`.
     */
    virtual void solve(Core::LinAlg::SparseMatrix& jacobian,
        const Core::LinAlg::Vector<double>& residual, const Core::LinAlg::Vector<double>& x,
        const NewtonLinearSystemMetadata& metadata, Core::LinAlg::Vector<double>& delta) = 0;
  };

  /**
   * @brief Context for constructing the baseline sparse Newton linear solver.
   */
  struct SparseNewtonLinearSolverContext
  {
    MPI_Comm comm;
    const Teuchos::ParameterList& linear_solver_parameters;
    std::function<const Teuchos::ParameterList&(int)> solver_params_callback;
    const Core::LinAlg::Map& correction_map;
    SparseNewtonLinearSolverProfile* profile = nullptr;
  };

  /**
   * @brief Reduced-lung Newton linear solver backed by Core::LinAlg::Solver.
   */
  class SparseNewtonLinearSolver : public NewtonLinearSolver
  {
   public:
    explicit SparseNewtonLinearSolver(const SparseNewtonLinearSolverContext& context);

    void solve(Core::LinAlg::SparseMatrix& jacobian, const Core::LinAlg::Vector<double>& residual,
        const Core::LinAlg::Vector<double>& x, const NewtonLinearSystemMetadata& metadata,
        Core::LinAlg::Vector<double>& delta) override;

   private:
    std::shared_ptr<Core::LinAlg::Solver> linear_solver_;
    std::unique_ptr<Core::LinAlg::Vector<double>> rhs_;
    SparseNewtonLinearSolverProfile* profile_ = nullptr;
  };
}  // namespace ReducedLung

FOUR_C_NAMESPACE_CLOSE

#endif
