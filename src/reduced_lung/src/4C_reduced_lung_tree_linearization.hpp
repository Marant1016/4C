// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_REDUCED_LUNG_TREE_LINEARIZATION_HPP
#define FOUR_C_REDUCED_LUNG_TREE_LINEARIZATION_HPP

#include "4C_config.hpp"

#include <span>
#include <utility>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace ReducedLung
{
  class TreeCoefficientAssemblyTarget
  {
   public:
    virtual ~TreeCoefficientAssemblyTarget() = default;

    virtual void append_value(int local_row_id, int local_dof_id, double value) = 0;

    virtual void replace_value(int local_row_id, int local_dof_id, double value) = 0;

    virtual void replace_values(std::span<const int> local_row_ids,
        std::span<const int> local_dof_ids, std::span<const double> values);
  };

  /**
   * @brief Structured local derivative storage for the tree Newton linear solver.
   *
   * Coefficients are addressed by local residual row id and local dof id from the reduced-lung
   * locally-relevant dof map. This mirrors the coefficient lookup that the validation tree path
   * currently performs on the completed sparse Jacobian, but avoids constructing that sparse
   * matrix for the production `NewtonTree` path.
   */
  class TreeLinearization : public TreeCoefficientAssemblyTarget
  {
   public:
    TreeLinearization() = default;

    TreeLinearization(int num_rows, int num_dofs) { reset(num_rows, num_dofs); }

    void reset(int num_rows, int num_dofs);

    void clear_values();

    void reserve_row_entries(int local_row_id, int entry_count);

    void set_value(int local_row_id, int local_dof_id, double value);

    void append_value(int local_row_id, int local_dof_id, double value) override;

    void replace_value(int local_row_id, int local_dof_id, double value) override;

    void replace_values(std::span<const int> local_row_ids, std::span<const int> local_dof_ids,
        std::span<const double> values) override;

    [[nodiscard]] double value(int local_row_id, int local_dof_id) const;

    [[nodiscard]] const std::vector<std::pair<int, double>>& entries(int local_row_id) const;

    [[nodiscard]] int num_rows() const { return num_rows_; }

    [[nodiscard]] int num_dofs() const { return num_dofs_; }

   private:
    int num_rows_ = 0;
    int num_dofs_ = 0;
    std::vector<std::vector<std::pair<int, double>>> rows_;
  };
}  // namespace ReducedLung

FOUR_C_NAMESPACE_CLOSE

#endif
