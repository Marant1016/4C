No code was modified. I ignored `src/reduced_lung/src/1d_pipe_flow/*` and traced the non-1D reduced-lung workflow.

**Big Picture**
The reduced lung model solves a nonlinear algebraic system each time step:

`F(x, time) = 0`

Here, `x` is the vector of reduced-lung unknowns, mainly pressures and flows. NOX drives the nonlinear Newton iteration. UMFPACK is only used inside each Newton step to solve the linearized system.

The split is:

| Level | Solver | Responsibility |
| --- | --- | --- |
| Nonlinear | NOX | Repeatedly evaluate residual `F(x)`, Jacobian `J(x)`, compute Newton updates, check convergence |
| Linear | UMFPACK via Amesos2 | Solve one sparse linear system per Newton step |

**Main Call Flow**
1. Entry point is `ReducedLung::reduced_lung_main()` in `src/reduced_lung/src/4C_reduced_lung_main.cpp:328`.

2. `make_reduced_lung_context_from_problem()` reads `reduced_dimensional_lung` parameters and the selected linear solver parameters from `Global::Problem` in `4C_reduced_lung_main.cpp:56`.

3. `ReducedLungSimulation::initialize()` builds the model in `4C_reduced_lung_main.cpp:93`.

4. `build_linear_system_and_solver()` creates vectors, sparse matrix, assembly callbacks, and `ReducedLung::NoxSolver` in `4C_reduced_lung_main.cpp:222`.

5. Each time step calls `nox_solver_->solve(current_time_)` in `ReducedLungSimulation::solve_timestep()` at `4C_reduced_lung_main.cpp:253`.

6. `ReducedLung::NoxSolver::solve()` calls `adapter_->solve()` in `4C_reduced_lung_helpers.cpp:168`.

7. `NOX::Nln::Adapter::solve()` runs NOX and copies the converged NOX solution back into the bound reduced-lung solution vector in `src/solver_nonlin_nox/4C_solver_nonlin_nox_adapter.cpp:146`.

8. After NOX converges, reduced-lung state is synchronized and end-of-step history variables are advanced in `4C_reduced_lung_main.cpp:271`.

**Where Reduced Lung Calls Solver Infrastructure**
The main reduced-lung solver handoff is in `ReducedLung::NoxSolver`.

Defined in:

`src/reduced_lung/src/4C_reduced_lung_helpers.hpp:144`

Implemented in:

`src/reduced_lung/src/4C_reduced_lung_helpers.cpp:123`

Important calls:

`4C_reduced_lung_main.cpp:250` constructs `NoxSolver`.

`4C_reduced_lung_main.cpp:269` calls `nox_solver_->solve(current_time_)`.

`4C_reduced_lung_helpers.cpp:130` constructs `Core::LinAlg::Solver`.

`4C_reduced_lung_helpers.cpp:164` constructs `NOX::Nln::Adapter`.

`4C_reduced_lung_helpers.cpp:171` calls `adapter_->solve()`.

So reduced lung itself does not manually run Newton iterations. It gives NOX the current solution vector, a Jacobian operator, a linear solver, and callbacks for residual/Jacobian assembly.

**NOX Role**
NOX is the nonlinear solver driver.

Reduced lung configures NOX in `create_nox_parameter_list()` at `4C_reduced_lung_helpers.cpp:451`.

Current settings:

| NOX setting | Value |
| --- | --- |
| Nonlinear solver | `"Line Search Based"` |
| Direction method | Newton |
| Line search | `"Full Step"` |
| Residual tolerance | `dynamics.nonlinear_residual_tolerance` |
| Increment tolerance | `dynamics.nonlinear_increment_tolerance` |
| Max nonlinear iterations | `dynamics.max_nonlinear_iterations` |

The Newton direction is handled by `NOX::Nln::Direction::Newton` in `src/solver_nonlin_nox/4C_solver_nonlin_nox_direction_newton.cpp:29`.

That direction first computes residual and Jacobian together:

`NOX::Nln::Direction::Newton::compute()` calls `nlnSoln->compute_f_and_jacobian()` at `4C_solver_nonlin_nox_direction_newton.cpp:42`.

Then the standard NOX Newton machinery solves for the Newton direction.

**NOX Inputs**
NOX receives these from reduced lung:

| Input | Source |
| --- | --- |
| Initial solution vector `x` | `context.x`, created as `x_` in `4C_reduced_lung_main.cpp:231` |
| Jacobian operator | `context.jacobian`, created as `sysmat_` in `4C_reduced_lung_main.cpp:232` |
| Residual callback | Bound to `NoxSolver::residual()` in `4C_reduced_lung_helpers.cpp:150` |
| Jacobian callback | Bound to `NoxSolver::jacobian()` in `4C_reduced_lung_helpers.cpp:154` |
| Linear solver map | `sol_generic -> linear_solver_` in `4C_reduced_lung_helpers.cpp:159` |
| NOX parameters | `create_nox_parameter_list()` in `4C_reduced_lung_helpers.cpp:451` |
| MPI communicator | `context.comm` |
| Current time | stored in `NoxSolver::current_time_` before each solve |

The first time step starts from the zero vector. Later time steps reuse the previous converged solution because `Adapter::solve()` copies the final NOX solution back into the bound vector at `4C_solver_nonlin_nox_adapter.cpp:164`.

**NOX Outputs**
NOX outputs:

| Output | Where |
| --- | --- |
| Converged solution vector | Copied back into `x_solution_` through `x_nox_` in `4C_solver_nonlin_nox_adapter.cpp:164` |
| Number of nonlinear iterations | Returned by `Adapter::solve()` at `4C_solver_nonlin_nox_adapter.cpp:165` |
| Failure if not converged | `NOX::Nln::Problem::check_final_status()` throws in `4C_solver_nonlin_nox_problem.cpp:161` |

Reduced lung then copies that converged `x` into its own dof vectors in `NoxSolver::sync_state_from_x()` at `4C_reduced_lung_helpers.cpp:211`.

**UMFPACK Role**
UMFPACK is the direct sparse linear solver used inside Newton.

Reduced lung selects the linear solver through:

`dynamics.linear_solver` in `src/reduced_lung/src/4C_reduced_lung_input.hpp:43`

The solver parameter list is obtained here:

`problem.solver_params(parameters.dynamics.linear_solver)` in `4C_reduced_lung_main.cpp:73`

`Core::LinAlg::Solver` translates that parameter list in `src/core/linear_solver/src/method/4C_linear_solver_method_linalg.cpp:399`.

If the solver type is `UMFPACK`, `Core::LinAlg::Solver::setup()` creates a `Core::LinearSolver::DirectSolver` in `4C_linear_solver_method_linalg.cpp:178`.

The direct solver maps UMFPACK to Amesos2 solver type `"Umfpack"` in:

`src/core/linear_solver/src/method/4C_linear_solver_method_direct.cpp:89`

UMFPACK does not know anything about lung physics. It only receives a sparse matrix and a right-hand side vector from NOX’s linear system layer.

**UMFPACK Inputs**
During each Newton step, NOX calls:

`NOX::Nln::LinearSystem::apply_jacobian_inverse()` in `src/solver_nonlin_nox/4C_solver_nonlin_nox_linearsystem.cpp:200`.

That constructs the linear problem:

`J * y = F`

Specifically:

| UMFPACK input | Source |
| --- | --- |
| Sparse matrix `J` | Reduced-lung Jacobian `sysmat_` |
| Right-hand side `F` | Current NOX residual vector |
| Output vector `y` | NOX Newton work vector before sign flip |

The call into the 4C linear solver is:

`currSolver->solve(linProblem.jac, linProblem.lhs, linProblem.rhs, solver_params)` at `4C_solver_nonlin_nox_linearsystem.cpp:271`.

Then `Core::LinAlg::Solver::solve()` calls the direct solver in `4C_linear_solver_method_linalg.cpp:222`.

Finally, `Core::LinearSolver::DirectSolver::solve()` calls Amesos2/UMFPACK:

`symbolicFactorization()`, `numericFactorization()`, then `solve()` in `4C_linear_solver_method_direct.cpp:146`.

**UMFPACK Outputs**
UMFPACK writes the solution of the linear system into the NOX result vector.

Important sign detail:

`apply_jacobian_inverse()` solves:

`J * y = F`

Then `GroupBase::computeNewton()` scales the result by `-1`:

`NewtonVector.scale(-1.0)` in `src/solver_nonlin_nox/4C_solver_nonlin_nox_group_base.cpp:194`.

So the actual Newton direction is:

`d = -J^{-1} F`

**Residual Assembly**
Residual assembly enters reduced lung through:

`ReducedLung::NoxSolver::residual()` in `4C_reduced_lung_helpers.cpp:177`.

The steps are:

1. `sync_state_from_x(x)` copies the NOX solution vector into reduced-lung dof vectors.

2. `Core::LinAlg::export_to(x, dofs_)` copies into locally owned dofs in `4C_reduced_lung_helpers.cpp:213`.

3. `Core::LinAlg::export_to(dofs_, locally_relevant_dofs_)` fills the ghosted/locally relevant vector in `4C_reduced_lung_helpers.cpp:214`.

4. State updaters run in `4C_reduced_lung_helpers.cpp:216`.

5. Residual assembler callbacks run in order in `4C_reduced_lung_helpers.cpp:182`.

The residual callback order is created by `create_default_nox_assembly_pipeline()` in `4C_reduced_lung_helpers.cpp:42`.

Order:

| Order | Component | Function |
| --- | --- | --- |
| 1 | Airways | `Airways::update_residual_vector()` |
| 2 | Terminal units | `TerminalUnits::update_residual_vector()` |
| 3 | Junctions | `Junctions::update_residual_vector()` |
| 4 | Boundary conditions | `BoundaryConditions::update_residual_vector()` |

The assemblers use `replace_local_value()`, meaning each equation row is directly overwritten with its current residual value.

Examples:

Airway residual rows are written in `src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.cpp:21`.

Terminal-unit residual rows are written in `src/reduced_lung/src/terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp:25`.

Junction residual rows are written in `src/reduced_lung/src/4C_reduced_lung_junctions.cpp:258`.

Boundary-condition residual rows are written in `src/reduced_lung/src/4C_reduced_lung_boundary_conditions.cpp:386`.

Conceptually, residual entries are:

| Component | Meaning |
| --- | --- |
| Airway rows | Element pressure-flow equations |
| Terminal-unit rows | Terminal-unit pressure/flow/elastic/rheology equations |
| Junction rows | Pressure continuity and mass conservation at connections/bifurcations |
| Boundary rows | Imposed pressure or flow conditions |

**Jacobian Assembly**
Jacobian assembly enters reduced lung through:

`ReducedLung::NoxSolver::jacobian()` in `4C_reduced_lung_helpers.cpp:190`.

The steps are:

1. `sync_state_from_x(x)` updates dofs and model internal state.

2. The generic NOX sparse operator is cast to `Core::LinAlg::SparseMatrix`.

3. Jacobian assembler callbacks run in order.

4. If the matrix is not filled yet, `jac_matrix->complete()` is called at `4C_reduced_lung_helpers.cpp:206`.

The Jacobian callback order is also set by `create_default_nox_assembly_pipeline()`:

| Order | Component | Function |
| --- | --- | --- |
| 1 | Airways | `Airways::update_jacobian()` |
| 2 | Terminal units | `TerminalUnits::update_jacobian()` |
| 3 | Junctions | `Junctions::update_jacobian()` |
| 4 | Boundary conditions | `BoundaryConditions::update_jacobian()` |

Important implementation detail:

On the first Jacobian call, most components insert sparse matrix entries using `insert_my_values()`.

On later Jacobian calls, the matrix is already filled. Components with nonlinear derivative values update existing entries using `replace_my_values()`.

Junction and boundary Jacobian entries are constant, so they return early if `sysmat.filled()` is true:

`Junctions::update_jacobian()` checks this at `4C_reduced_lung_junctions.cpp:306`.

`BoundaryConditions::update_jacobian()` checks this at `4C_reduced_lung_boundary_conditions.cpp:399`.

**Solution Update**
There are two update levels.

At the linear level:

`UMFPACK` computes `y` from:

`J * y = F`

NOX then flips the sign:

`d = -y`

This is done in `GroupBase::computeNewton()` at `4C_solver_nonlin_nox_group_base.cpp:180`.

At the nonlinear level:

NOX applies the Newton direction to the current nonlinear solution.

The update is performed by `NOX::Nln::Group::computeX()`:

`x_new = x_old + step * d`

This is in `src/solver_nonlin_nox/4C_solver_nonlin_nox_group.cpp:85`.

Because reduced lung configures `"Full Step"` line search in `4C_reduced_lung_helpers.cpp:460`, the normal step size is `1.0`.

After NOX convergence:

`NOX::Nln::Adapter::solve()` copies the final NOX group solution into the bound reduced-lung `x` vector in `4C_solver_nonlin_nox_adapter.cpp:164`.

`ReducedLung::NoxSolver::solve()` then calls `sync_state_from_x(x_solution_)` in `4C_reduced_lung_helpers.cpp:172`.

Then the time-step routine advances history variables:

`TerminalUnits::end_of_timestep_routine()` in `4C_reduced_lung_main.cpp:271`.

`Airways::end_of_timestep_routine()` in `4C_reduced_lung_main.cpp:272`.

**Important Classes And Files**
| Area | Main files/classes/functions |
| --- | --- |
| Simulation driver | `4C_reduced_lung_main.cpp`, `ReducedLungSimulation`, `initialize()`, `run()`, `solve_timestep()` |
| Reduced-lung NOX wrapper | `4C_reduced_lung_helpers.hpp/cpp`, `NoxSolver`, `NoxSolverContext`, `NoxAssemblyPipeline` |
| NOX adapter | `4C_solver_nonlin_nox_adapter.hpp/cpp`, `NOX::Nln::Adapter` |
| NOX problem/group | `4C_solver_nonlin_nox_problem.cpp`, `4C_solver_nonlin_nox_group.cpp`, `4C_solver_nonlin_nox_group_base.cpp` |
| NOX Newton direction | `4C_solver_nonlin_nox_direction_newton.cpp` |
| NOX linear system | `4C_solver_nonlin_nox_linearsystem.cpp`, `4C_solver_nonlin_nox_linearsystem_generic.cpp` |
| 4C linear solver wrapper | `4C_linear_solver_method_linalg.hpp/cpp`, `Core::LinAlg::Solver` |
| UMFPACK direct solver | `4C_linear_solver_method_direct.hpp/cpp`, `Core::LinearSolver::DirectSolver` |
| Airway assembly | `airways/4C_reduced_lung_airways.cpp`, `airways/4C_reduced_lung_airways_wall_mechanics.cpp`, `airways/4C_reduced_lung_airways_flow_resistance.cpp` |
| Terminal-unit assembly | `terminal_units/4C_reduced_lung_terminal_unit.cpp`, `terminal_units/4C_reduced_lung_terminal_unit_rheology.cpp`, `terminal_units/4C_reduced_lung_terminal_unit_elasticity.cpp` |
| Junction assembly | `4C_reduced_lung_junctions.hpp/cpp` |
| Boundary-condition assembly | `4C_reduced_lung_boundary_conditions.hpp/cpp` |

**Where A Tree-Based Solver Could Fit Later**
The cleanest insertion point depends on what the tree solver replaces.

If the new tree solver replaces the whole nonlinear solve, insert it beside or instead of `NoxSolver` in:

`ReducedLungSimulation::build_linear_system_and_solver()` at `4C_reduced_lung_main.cpp:222`

and call it from:

`ReducedLungSimulation::solve_timestep()` at `4C_reduced_lung_main.cpp:253`.

If the new tree solver is only a linear solver for Newton corrections, keep NOX and replace the linear backend below:

`Core::LinAlg::Solver` in `4C_linear_solver_method_linalg.cpp`

or the NOX linear system layer in:

`4C_solver_nonlin_nox_linearsystem.cpp`

If the new tree solver still wants the existing residual/Jacobian assembly, reuse:

`NoxAssemblyPipeline` in `4C_reduced_lung_helpers.hpp:67`

and possibly create a new solver class with the same high-level interface as `NoxSolver::solve(time)`.

Most likely practical path:

Create a new reduced-lung solver class parallel to `NoxSolver`, selected in `build_linear_system_and_solver()`, that reuses the same model containers, maps, dof vectors, and end-of-timestep routines.

**Summary Diagram**
```text
ReducedLungSimulation::solve_timestep()
  -> NoxSolver::solve(time)
    -> NOX::Nln::Adapter::solve()
      -> NOX nonlinear iteration
        -> residual callback: NoxSolver::residual(x, F)
          -> sync x into dof vectors
          -> assemble airway + terminal + junction + BC residuals
        -> jacobian callback: NoxSolver::jacobian(x, J)
          -> sync x into dof vectors
          -> assemble airway + terminal + junction + BC Jacobian
        -> NOX linear system
          -> Core::LinAlg::Solver::solve(J, y, F)
            -> DirectSolver(UMFPACK)
            -> solve J y = F
        -> NOX forms Newton direction d = -y
        -> NOX updates x_new = x_old + d
      -> repeat until residual/update tolerance or failure
    -> copy converged NOX x back to reduced-lung x
    -> sync reduced-lung dof/state vectors
  -> terminal-unit and airway end-of-timestep history updates
```
