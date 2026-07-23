// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_reduced_lung_newton_solver.hpp"

#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_reduced_lung_solver_profile.hpp"
#include "4C_reduced_lung_tree_linearization.hpp"
#include "4C_utils_exceptions.hpp"

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

  NewtonSolver::NewtonSolver(const NewtonSolverContext& context, double initial_time)
      : x_solution_(context.x),
        dofs_(context.dofs),
        locally_relevant_dofs_(context.locally_relevant_dofs),
        jacobian_(context.jacobian),
        assembly_pipeline_(context.assembly_pipeline),
        residual_(context.x.get_map(), true),
        delta_(context.x.get_map(), true),
        tree_linearization_(context.x.local_length(), context.locally_relevant_dofs.local_length()),
        dt_(context.dynamics.time_increment),
        current_time_(initial_time),
        max_nonlinear_iterations_(
            static_cast<unsigned int>(context.dynamics.max_nonlinear_iterations)),
        nonlinear_residual_tolerance_(context.dynamics.nonlinear_residual_tolerance),
        nonlinear_increment_tolerance_(context.dynamics.nonlinear_increment_tolerance),
        linear_solver_(context.linear_solver),
        profile_(context.profile)
  {
    if (context.dynamics.max_nonlinear_iterations <= 0)
    {
      FOUR_C_THROW(
          "ReducedLung::NewtonSolver requires a positive max_nonlinear_iterations, got {}.",
          context.dynamics.max_nonlinear_iterations);
    }
    if (linear_solver_ == nullptr)
    {
      FOUR_C_THROW("ReducedLung::NewtonSolver requires a valid Newton linear solver instance.");
    }
    if (assembly_pipeline_.residual_assemblers.empty())
    {
      FOUR_C_THROW("ReducedLung::NewtonSolver requires at least one residual assembler callback.");
    }
    if (assembly_pipeline_.jacobian_assemblers.empty())
    {
      FOUR_C_THROW("ReducedLung::NewtonSolver requires at least one Jacobian assembler callback.");
    }
    if (linear_solver_->linearization_type() == NewtonLinearizationType::StructuredTreeBlocks &&
        assembly_pipeline_.tree_linearization_assemblers.empty())
    {
      FOUR_C_THROW(
          "ReducedLung::NewtonSolver requires tree-linearization assemblers for structured tree "
          "linear solves.");
    }
  }

  unsigned int NewtonSolver::solve(double time)
  {
    const auto solve_start = Clock::now();
    current_time_ = time;
    double increment_norm = 0.0;
    if (profile_ != nullptr)
    {
      profile_->last_residual_norms.clear();
      profile_->last_increment_norms.clear();
    }

    for (unsigned int iteration = 0; iteration <= max_nonlinear_iterations_; ++iteration)
    {
      const auto sync_start = Clock::now();
      sync_state_from_x(x_solution_);
      if (profile_ != nullptr)
      {
        profile_->state_sync_time += elapsed_seconds(sync_start);
      }
      const double residual_norm = assemble_residual_for_current_state();
      if (profile_ != nullptr)
      {
        profile_->last_residual_norms.push_back(residual_norm);
      }
      const bool residual_converged = residual_norm <= nonlinear_residual_tolerance_;
      const bool increment_converged =
          iteration == 0 || increment_norm <= nonlinear_increment_tolerance_;

      if (residual_converged && increment_converged)
      {
        if (profile_ != nullptr)
        {
          profile_->total_solve_time += elapsed_seconds(solve_start);
          profile_->last_nonlinear_iterations = iteration;
          ++profile_->solve_count;
        }
        return iteration;
      }

      if (iteration == max_nonlinear_iterations_)
      {
        FOUR_C_THROW(
            "ReducedLung::NewtonSolver did not converge at time {} after {} Newton corrections. "
            "Final residual norm: {}, final increment norm: {}.",
            current_time_, max_nonlinear_iterations_, residual_norm, increment_norm);
      }

      if (linear_solver_->linearization_type() == NewtonLinearizationType::StructuredTreeBlocks)
      {
        assemble_tree_linearization_for_current_state();
      }
      else
      {
        assemble_jacobian_for_current_state();
      }
      increment_norm = solve_linear_correction(iteration);
      if (profile_ != nullptr)
      {
        profile_->last_increment_norms.push_back(increment_norm);
      }
      x_solution_.update(1.0, delta_, 1.0);
    }

    FOUR_C_THROW("ReducedLung::NewtonSolver reached an unreachable nonlinear-solver state.");
  }

  void NewtonSolver::sync_state_from_x(const Core::LinAlg::Vector<double>& x)
  {
    Core::LinAlg::export_to(x, dofs_);
    Core::LinAlg::export_to(dofs_, locally_relevant_dofs_);

    for (const auto& update_state : assembly_pipeline_.state_updaters)
    {
      update_state(locally_relevant_dofs_, dt_);
    }
  }

  double NewtonSolver::assemble_residual_for_current_state()
  {
    residual_.put_scalar(0.0);
    const auto assembly_start = Clock::now();
    for (const auto& assemble_residual : assembly_pipeline_.residual_assemblers)
    {
      assemble_residual(residual_, locally_relevant_dofs_, current_time_, dt_);
    }
    if (profile_ != nullptr)
    {
      profile_->residual_assembly_time += elapsed_seconds(assembly_start);
    }

    double residual_norm = 0.0;
    residual_.norm_2(&residual_norm);
    return residual_norm;
  }

  void NewtonSolver::assemble_jacobian_for_current_state()
  {
    const auto assembly_start = Clock::now();
    for (const auto& assemble_jacobian : assembly_pipeline_.jacobian_assemblers)
    {
      assemble_jacobian(jacobian_, locally_relevant_dofs_, current_time_, dt_);
    }
    if (profile_ != nullptr)
    {
      profile_->sparse_jacobian_assembly_time += elapsed_seconds(assembly_start);
    }

    if (!jacobian_.filled())
    {
      const auto complete_start = Clock::now();
      jacobian_.complete();
      if (profile_ != nullptr)
      {
        profile_->sparse_jacobian_complete_time += elapsed_seconds(complete_start);
      }
    }
  }

  void NewtonSolver::assemble_tree_linearization_for_current_state()
  {
    const auto assembly_start = Clock::now();
    tree_linearization_.reset(residual_.local_length(), locally_relevant_dofs_.local_length());
    for (const auto& assemble_tree_linearization : assembly_pipeline_.tree_linearization_assemblers)
    {
      assemble_tree_linearization(tree_linearization_, locally_relevant_dofs_, current_time_, dt_);
    }
    linear_solver_->set_tree_linearization(tree_linearization_);
    if (profile_ != nullptr)
    {
      profile_->structured_tree_linearization_assembly_time += elapsed_seconds(assembly_start);
    }
  }

  double NewtonSolver::solve_linear_correction(unsigned int iteration)
  {
    const NewtonLinearSystemMetadata metadata{
        .current_time = current_time_,
        .time_step_size_dt = dt_,
        .nonlinear_iteration = iteration,
    };
    const auto linear_solve_start = Clock::now();
    linear_solver_->solve(jacobian_, residual_, x_solution_, metadata, delta_);
    if (profile_ != nullptr)
    {
      profile_->linear_solve_time += elapsed_seconds(linear_solve_start);
    }

    double increment_norm = 0.0;
    delta_.norm_2(&increment_norm);
    return increment_norm;
  }
}  // namespace ReducedLung

FOUR_C_NAMESPACE_CLOSE
