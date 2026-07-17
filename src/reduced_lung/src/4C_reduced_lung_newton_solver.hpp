// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_REDUCED_LUNG_NEWTON_SOLVER_HPP
#define FOUR_C_REDUCED_LUNG_NEWTON_SOLVER_HPP

#include "4C_config.hpp"

#include "4C_reduced_lung_helpers.hpp"

#include <mpi.h>
#include <Teuchos_ParameterList.hpp>

#include <functional>
#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class Solver;
  class SparseMatrix;
  template <typename T>
  class Vector;
}  // namespace Core::LinAlg

namespace ReducedLung
{
  /**
   * @brief Context bundling all objects required by @ref NewtonSolver.
   */
  struct NewtonSolverContext
  {
    MPI_Comm comm;  ///< MPI communicator used by the linear solver.
    const ReducedLungParameters::Dynamics& dynamics;  ///< Nonlinear/timestep solver parameters.
    const Teuchos::ParameterList& linear_solver_parameters;  ///< Linear solver configuration.
    std::function<const Teuchos::ParameterList&(int)>
        solver_params_callback;  ///< Callback for nested/ID-based solver parameters.
    const ReducedLungAssemblyPipeline& assembly_pipeline;  ///< Ordered model assembly callbacks.
    Core::LinAlg::Vector<double>& dofs;                    ///< Owned dof vector.
    Core::LinAlg::Vector<double>& locally_relevant_dofs;   ///< Ghosted dof vector.
    Core::LinAlg::Vector<double>& x;                       ///< Nonlinear solution vector.
    Core::LinAlg::SparseMatrix& jacobian;                  ///< Newton-system Jacobian matrix.
  };

  /**
   * @brief Full-step Newton solver for reduced-lung nonlinear systems.
   *
   * This solver is intentionally parallel to @ref NoxSolver. It reuses the reduced-lung assembly
   * pipeline and the existing sparse linear solver backend, but owns the nonlinear iteration loop.
   */
  class NewtonSolver
  {
   public:
    /**
     * @brief Construct a custom Newton solver for reduced-lung systems.
     *
     * @param context Solver setup context including communicator, dynamics, linear solver setup,
     * assembly pipeline callbacks, and all bound vectors/matrices.
     * @param initial_time Initial time for the simulation.
     */
    NewtonSolver(const NewtonSolverContext& context, double initial_time = 0.0);

    NewtonSolver(const NewtonSolver&) = delete;
    NewtonSolver& operator=(const NewtonSolver&) = delete;
    NewtonSolver(NewtonSolver&&) = delete;
    NewtonSolver& operator=(NewtonSolver&&) = delete;
    ~NewtonSolver() = default;

    /**
     * @brief Solve the nonlinear system at the given physical time.
     *
     * @param time Current physical time for time-dependent boundary conditions.
     * @return Number of Newton corrections applied.
     */
    unsigned int solve(double time);

   private:
    void sync_state_from_x(const Core::LinAlg::Vector<double>& x);

    double assemble_residual_for_current_state();

    void assemble_jacobian_for_current_state();

    double solve_linear_correction(unsigned int iteration);

    Core::LinAlg::Vector<double>& x_solution_;
    Core::LinAlg::Vector<double>& dofs_;
    Core::LinAlg::Vector<double>& locally_relevant_dofs_;
    Core::LinAlg::SparseMatrix& jacobian_;
    ReducedLungAssemblyPipeline assembly_pipeline_;

    Core::LinAlg::Vector<double> residual_;
    Core::LinAlg::Vector<double> rhs_;
    Core::LinAlg::Vector<double> delta_;

    double dt_;
    double current_time_;
    unsigned int max_nonlinear_iterations_;
    double nonlinear_residual_tolerance_;
    double nonlinear_increment_tolerance_;

    std::shared_ptr<Core::LinAlg::Solver> linear_solver_;
  };
}  // namespace ReducedLung

FOUR_C_NAMESPACE_CLOSE

#endif
