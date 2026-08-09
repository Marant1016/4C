// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_REDUCED_LUNG_TREE_METADATA_HPP
#define FOUR_C_REDUCED_LUNG_TREE_METADATA_HPP

#include "4C_config.hpp"

#include "4C_reduced_lung_airways.hpp"
#include "4C_reduced_lung_boundary_conditions.hpp"
#include "4C_reduced_lung_input.hpp"
#include "4C_reduced_lung_junctions.hpp"
#include "4C_reduced_lung_terminal_unit.hpp"

#include <array>
#include <cstdint>
#include <map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class Map;
}

namespace ReducedLung
{
  enum class TreeElementKind : std::uint8_t
  {
    Airway,
    TerminalUnit
  };

  enum class TreeJunctionKind : std::uint8_t
  {
    Connection,
    Bifurcation
  };

  enum class TreeBoundarySide : std::uint8_t
  {
    Inlet,
    Outlet
  };

  /**
   * @brief Per-element topology, dof, and state-equation layout for tree solvers.
   */
  struct TreeElementMetadata
  {
    int global_element_id = -1;
    TreeElementKind kind = TreeElementKind::Airway;

    int inlet_node_id = -1;
    int outlet_node_id = -1;

    int parent_element_index = -1;
    std::array<int, 2> child_element_indices{-1, -1};
    int child_count = 0;

    int first_global_dof = -1;
    int num_dofs = 0;
    std::vector<int> global_dof_ids;
    std::vector<int> local_dof_ids;

    int first_local_state_equation_id = -1;
    int first_global_state_equation_id = -1;
    int num_state_equations = 0;
    int owner_rank = -1;

    [[nodiscard]] bool is_leaf() const { return child_count == 0; }
  };

  /**
   * @brief Junction-equation metadata connecting one parent element to one or two children.
   */
  struct TreeJunctionMetadata
  {
    TreeJunctionKind kind = TreeJunctionKind::Connection;
    int parent_element_index = -1;
    std::array<int, 2> child_element_indices{-1, -1};
    int child_count = 0;

    int first_local_equation_id = -1;
    int first_global_equation_id = -1;
    int num_equations = 0;
    int owner_rank = -1;

    std::vector<int> global_dof_ids;
    std::vector<int> local_dof_ids;
  };

  /**
   * @brief Boundary-condition equation metadata attached to an inlet or outlet element side.
   */
  struct TreeBoundaryConditionMetadata
  {
    BoundaryConditions::Type type = BoundaryConditions::Type::Pressure;
    TreeBoundarySide side = TreeBoundarySide::Inlet;

    int node_id = -1;
    int element_index = -1;

    int local_equation_id = -1;
    int global_equation_id = -1;

    int global_dof_id = -1;
    int local_dof_id = -1;
    int owner_rank = -1;
  };

  /**
   * @brief Immutable tree metadata consumed by future tree-based Newton linear solvers.
   */
  struct ReducedLungTreeMetadata
  {
    std::vector<TreeElementMetadata> elements;
    std::vector<TreeJunctionMetadata> junctions;
    std::vector<TreeBoundaryConditionMetadata> boundary_conditions;

    std::map<int, int> element_index_by_global_id;

    int root_element_index = -1;
    int root_node_id = -1;

    std::vector<int> airway_element_indices;
    std::vector<int> terminal_unit_element_indices;

    std::vector<std::vector<int>> bottom_up_layers;
    std::vector<std::vector<int>> top_down_layers;

    int num_global_dofs = 0;
    int num_global_equations = 0;
    int num_locally_relevant_dofs = 0;
  };

  /**
   * @brief Inputs required to build reduced-lung tree metadata from the existing model setup.
   */
  struct ReducedLungTreeMetadataContext
  {
    const ReducedLungParameters& parameters;
    const std::map<int, int>& first_global_dof_of_ele;
    const std::map<int, int>& global_dof_per_ele;

    const Airways::AirwayContainer& airways;
    const TerminalUnits::TerminalUnitContainer& terminal_units;

    const Junctions::ConnectionData& connections;
    const Junctions::BifurcationData& bifurcations;

    const BoundaryConditions::BoundaryConditionContainer& boundary_conditions;

    const Core::LinAlg::Map& row_map;
    const Core::LinAlg::Map& locally_relevant_dof_map;
  };

  /**
   * @brief Build and validate tree metadata without changing solver state or assembly behavior.
   */
  ReducedLungTreeMetadata build_reduced_lung_tree_metadata(
      const ReducedLungTreeMetadataContext& context);
}  // namespace ReducedLung

FOUR_C_NAMESPACE_CLOSE

#endif
