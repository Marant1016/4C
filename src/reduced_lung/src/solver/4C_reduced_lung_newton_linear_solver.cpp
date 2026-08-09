// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_reduced_lung_newton_linear_solver.hpp"

#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_reduced_lung_solver_profiles.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_shared_ptr_from_ref.hpp"

#include <chrono>

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
  }  // namespace

  SparseNewtonLinearSolver::SparseNewtonLinearSolver(const SparseNewtonLinearSolverContext& context)
      : linear_solver_(std::make_shared<Core::LinAlg::Solver>(context.linear_solver_parameters,
            context.comm, context.solver_params_callback, Core::IO::Verbositylevel::minimal)),
        rhs_(std::make_unique<Core::LinAlg::Vector<double>>(context.correction_map, true)),
        profile_(context.profile)
  {
    if (!linear_solver_)
    {
      FOUR_C_THROW(
          "ReducedLung::SparseNewtonLinearSolver requires a valid linear solver instance.");
    }
  }

  void SparseNewtonLinearSolver::solve(Core::LinAlg::SparseMatrix& jacobian,
      const Core::LinAlg::Vector<double>& residual, const Core::LinAlg::Vector<double>& x,
      const NewtonLinearSystemMetadata& metadata, Core::LinAlg::Vector<double>& delta)
  {
    (void)x;

    const auto solve_start = Clock::now();
    rhs_->scale(-1.0, residual);
    delta.put_scalar(0.0);

    Core::LinAlg::SolverParams solver_params;
    solver_params.refactor = true;
    solver_params.reset = metadata.nonlinear_iteration == 0;
    if (linear_solver_->params().isParameter("Projector"))
    {
      solver_params.projector =
          linear_solver_->params().get<std::shared_ptr<Core::LinAlg::LinearSystemProjector>>(
              "Projector");
    }

    const int linear_solver_status = linear_solver_->solve(
        Core::Utils::shared_ptr_from_ref(jacobian), Core::Utils::shared_ptr_from_ref(delta),
        Core::Utils::shared_ptr_from_ref(*rhs_), solver_params);
    if (linear_solver_status != 0)
    {
      FOUR_C_THROW(
          "ReducedLung::SparseNewtonLinearSolver failed at time {}, Newton iteration {} "
          "with status {}.",
          metadata.current_time, metadata.nonlinear_iteration, linear_solver_status);
    }
    if (profile_ != nullptr)
    {
      profile_->solve_time += elapsed_seconds(solve_start);
      ++profile_->solve_count;
    }
  }
}  // namespace ReducedLung

FOUR_C_NAMESPACE_CLOSE
