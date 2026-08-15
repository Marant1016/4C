# reduced_lung Tree Solver Instructions

## Scope

- These instructions apply to `src/reduced_lung/src/solver/tree/` and supplement the parent reduced-lung and solver instructions.
- This directory owns metadata, structured coefficient storage, and the serial tree-based Newton correction linear solver used by `NewtonTree`.

## Tree Solver Context

- The tree solver solves the same Newton correction system as the sparse path: `J * delta = -F`, followed by the nonlinear solver applying `x += delta`.
- The optimized algorithm condenses subtrees bottom-up into affine inlet relations `delta_q_in = G * delta_p_in + h`, closes the root, then recovers corrections top-down.
- `NewtonTree` is serial-only. Do not add silent fallback to `Nox` or `NewtonSparse` for multi-rank runs because that hides invalid performance results.
- Runtime `NewtonTree` should use `StructuredTreeBlocks`. `SparseJacobian` coefficient source is a validation/debugging path, not the normal production path.

## Important Files

- `4C_reduced_lung_tree_metadata.hpp` and `4C_reduced_lung_tree_metadata.cpp`: tree topology, dof/equation metadata, traversal layers, validation, and solver assumptions.
- `4C_reduced_lung_tree_linearization.hpp` and `4C_reduced_lung_tree_linearization.cpp`: row-oriented structured coefficient storage and `TreeCoefficientAssemblyTarget` interface.
- `4C_reduced_lung_tree_linear_solver.hpp` and `4C_reduced_lung_tree_linear_solver.cpp`: symbolic plan, direct coefficient target, bottom-up condensation, root closure, top-down recovery, batching, fallbacks, and profiling.
- `CMakeLists.txt`: this module uses `four_c_auto_define_module()`; keep new files aligned with the module auto-discovery pattern.

## Metadata Invariants

- Treat input topology as directional. In `element_nodes`, the first node is the inlet or parent side and the second node is the outlet or child side.
- Preserve validation for a single rooted directed tree, all elements reachable from the root, no cycles, branch degree at most two, terminal units as leaves, square system layout, a root inlet boundary, and one outlet boundary per leaf.
- Do not change dof ordering or row ordering from the reduced-lung module as part of tree-solver work.
- Keep connection and bifurcation metadata consistent with the topology. A metadata mismatch should fail fast rather than being repaired silently.

## Structured Coefficients

- Keep derivative ownership in the physics modules. Tree solver code should consume coefficients through `TreeCoefficientAssemblyTarget` or `TreeLinearization`.
- Static coefficients may be appended once and dynamic coefficients replaced each nonlinear iteration. Preserve this distinction for performance.
- Keep `TreeNewtonLinearSolver::direct_tree_coefficient_target()` as the fast path for runtime structured assembly unless the task explicitly changes that interface.
- When adding coefficient locations, ensure local row ids and locally relevant dof ids are used consistently.

## Numerical And Performance Rules

- Preserve robust pivot checks and dense fallback behavior unless the task is explicitly about changing fallback policy.
- Production block sizes are currently `2x2` for rigid airways and terminal units and `3x3` for Kelvin-Voigt airways. Keep optimized paths focused on these common shapes unless new physics requires more.
- Keep generic dense fallback available for unsupported or numerically difficult element blocks unless validation shows it is safe to remove.
- Avoid sparse-matrix access in the runtime structured path. Sparse lookup is acceptable for validation and targeted debugging.
- Keep profiling counters meaningful for large serial cases, especially solve time, dense fallback counts, unsupported fallback counts, SIMD/group counters, element counts, and max block size.

## Verification

- For numerical changes, compare tree corrections against sparse corrections on small cases before relying on large tree benchmarks.
- Useful builds from the source root include `cmake --build build/debug --target reduced_lung_objs --parallel 4` and reduced-lung unit test targets when available.
- Useful tests include `ctest -R "^unittests_reduced_lung(\.np2)?$" --output-on-failure` and registered `NewtonTree` reduced-lung input tests.
- For performance runs, use a release executable, keep `NewtonTree` serial, and enable `FOUR_C_REDUCED_LUNG_TREE_PROFILE=1` when profiling tree phases.
