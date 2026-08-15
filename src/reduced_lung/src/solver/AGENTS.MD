# reduced_lung Solver Instructions

## Scope

- These instructions apply to `src/reduced_lung/src/solver/` and supplement `src/reduced_lung/AGENTS.md` plus the repository-level guidance.
- This directory owns the custom reduced-lung Newton solver shared by `NewtonSparse` and `NewtonTree`.
- The `tree/` subdirectory has additional tree-solver-specific instructions in its own `AGENTS.md`.

## Solver Context

- `NewtonSolver` is a full-step nonlinear Newton loop that runs parallel to the existing NOX workflow.
- `SparseNewtonLinearSolver` is the baseline custom Newton linear solver and delegates to `Core::LinAlg::Solver` using the sparse Jacobian.
- `TreeNewtonLinearSolver` is selected through the same `NewtonLinearSolver` interface but consumes structured tree blocks instead of requiring runtime sparse Jacobian assembly.
- The custom Newton sign convention is `J * delta = -F`, followed by `x += delta`. Keep this convention consistent across sparse and tree paths.
- The residual norm is the authoritative convergence check. Do not accept a small increment as convergence if the residual is still above tolerance.

## Important Files

- `4C_reduced_lung_newton_solver.hpp` and `4C_reduced_lung_newton_solver.cpp`: nonlinear iteration, state synchronization, residual/Jacobian/tree-linearization assembly dispatch, correction solve, convergence, and profiling.
- `4C_reduced_lung_newton_linear_solver.hpp` and `4C_reduced_lung_newton_linear_solver.cpp`: linear-solver interface, linearization type selection, sparse backend, and direct tree coefficient target hook.
- `4C_reduced_lung_solver_profiles.hpp`: profile data shared by custom Newton, sparse Newton, and tree linear solver paths.
- `tree/`: serial optimized tree-based Newton correction solver and tree metadata/linearization support.
- `CMakeLists.txt`: this module uses `four_c_auto_define_module()`; keep file naming and include paths consistent with 4C module discovery.

## Assembly And State Rules

- Preserve the high-level order in each nonlinear iteration: synchronize `x`, update model state, assemble residual, check residual convergence, assemble the selected linearization, solve correction, then update `x`.
- `dofs`, `locally_relevant_dofs`, and `x` have different roles. Keep ownership and ghost synchronization explicit before changing vector updates.
- Runtime `NewtonSparse` should assemble the sparse Jacobian and call the sparse linear solver.
- Runtime `NewtonTree` should assemble structured tree coefficients and call the tree linear solver. Avoid adding sparse Jacobian assembly to the production tree path unless the task is specifically validation or debugging.
- Keep physics residuals and derivative formulas in the airway, terminal-unit, junction, and boundary-condition modules. Solver code should orchestrate assembly, not duplicate physics formulas.
- Be careful with sparse matrix fill state. Existing sparse Jacobian assemblers may insert on an unfilled matrix and replace or skip entries after `SparseMatrix::complete()`.

## Error Handling And Profiling

- Use `FOUR_C_THROW` for invalid solver state, unsupported configurations, convergence failure, and numerical breakdowns.
- Keep error messages specific enough to identify the nonlinear iteration, time, residual norm, increment norm, or solver phase when available.
- Profile instrumentation should remain optional and low overhead when the profile pointer is null.
- Preserve reduced-lung tree profiling behavior controlled by `FOUR_C_REDUCED_LUNG_TREE_PROFILE=1`.

## Verification

- Prefer targeted builds from the 4C source root, for example `cmake --build build/debug --target reduced_lung_objs --parallel 4`.
- For behavior changes, compare `NewtonSparse` against the existing NOX path or compare `NewtonTree` against `NewtonSparse` on small cases when practical.
- Useful targeted tests include `ctest -R "^unittests_reduced_lung(\.np2)?$" --output-on-failure` and registered reduced-lung input tests matching the changed solver path.
- If local dependencies or build directories are unavailable, state the verification gap clearly.
