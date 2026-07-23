// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_REDUCED_LUNG_SOLVER_PROFILE_HPP
#define FOUR_C_REDUCED_LUNG_SOLVER_PROFILE_HPP

#include "4C_config.hpp"

#include <cstdint>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace ReducedLung
{
  struct NoxSolverProfile
  {
    double total_solve_time = 0.0;
    double state_sync_time = 0.0;
    double residual_assembly_time = 0.0;
    double sparse_jacobian_assembly_time = 0.0;
    double sparse_jacobian_complete_time = 0.0;
    unsigned int solve_count = 0;
    unsigned int residual_evaluation_count = 0;
    unsigned int jacobian_evaluation_count = 0;
    unsigned int last_nonlinear_iterations = 0;
  };

  struct NewtonSolverProfile
  {
    double total_solve_time = 0.0;
    double state_sync_time = 0.0;
    double residual_assembly_time = 0.0;
    double sparse_jacobian_assembly_time = 0.0;
    double sparse_jacobian_complete_time = 0.0;
    double structured_tree_linearization_assembly_time = 0.0;
    double linear_solve_time = 0.0;
    unsigned int solve_count = 0;
    unsigned int last_nonlinear_iterations = 0;
    std::vector<double> last_residual_norms;
    std::vector<double> last_increment_norms;
  };

  struct SparseNewtonLinearSolverProfile
  {
    double solve_time = 0.0;
    unsigned int solve_count = 0;
  };

  struct TreeNewtonLinearSolverProfile
  {
    double total_solve_time = 0.0;
    double bottom_up_time = 0.0;
    double top_down_time = 0.0;
    double dense_solve_time = 0.0;
    double coefficient_lookup_time = 0.0;
    double communication_time = 0.0;
    unsigned int solve_count = 0;
    std::uint64_t dense_solve_count = 0;
    std::uint64_t coefficient_lookup_count = 0;
    std::uint64_t communicated_coefficient_count = 0;
    std::uint64_t communicated_residual_count = 0;
    std::uint64_t communication_bytes = 0;
    std::uint64_t element_count = 0;
    std::uint64_t total_local_block_dofs = 0;
    int max_local_block_size = 0;
  };
}  // namespace ReducedLung

FOUR_C_NAMESPACE_CLOSE

#endif
