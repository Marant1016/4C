# Assemply Step 4

## Goal

Avoid temporary heap allocations in airway model coefficient evaluators by changing flow-resistance and inertia callbacks from vector-returning APIs to output-buffer APIs.

This step keeps the same residual, Jacobian, and structured tree-linearization coefficient values. It only changes where temporary derivative/value storage lives.

## Motivation

The rigid airway tree-linearization path previously allocated temporary vectors every time it evaluated resistance and inertia derivatives:

```cpp
auto resistance_derivative =
    resistance_derivative_evaluator(airway_data, locally_relevant_dofs, dt);

auto inertia_derivative = ...;
```

For the gen16 rigid case, this happens for many airway elements across many Newton tree assemblies. The vector allocations are avoidable because each evaluator already knows the required output size.

## Code Changes

### Flow Evaluator APIs Now Fill Spans

Files:

```text
src/reduced_lung/src/airways/4C_reduced_lung_airways_flow_resistance.hpp
src/reduced_lung/src/airways/4C_reduced_lung_airways_flow_resistance.cpp
```

Changed evaluator callback types from returning `std::vector<double>` or `std::pair<std::vector<double>, std::vector<double>>` to writing into caller-provided `std::span<double>` buffers.

Examples:

```cpp
using InertiaEvaluator = std::function<void(
    const AirwayData&, const std::vector<double>& area, std::span<double> inertia)>;
```

```cpp
using FlowResistanceDerivativeEvaluatorRigid = std::function<void(
    const AirwayData&, const Core::LinAlg::Vector<double>&, double,
    std::span<double>)>;
```

Kelvin-Voigt two-component derivative evaluators now fill two output spans:

```cpp
std::span<double> q1_derivative
std::span<double> q2_derivative
```

### Removed Internal Poiseuille Temporary Vectors

`ComputePoiseuilleResistance` now fills a provided output span instead of returning a newly allocated vector.

Nonlinear derivative paths compute per-element Poiseuille resistance directly inside their loops where only one scalar is needed.

### Reusable Scratch In Wall Mechanics Evaluators

File:

```text
src/reduced_lung/src/airways/4C_reduced_lung_airways_wall_mechanics.cpp
```

The generated residual, Jacobian, and structured tree-linearization evaluator lambdas now capture reusable scratch vectors:

```cpp
return [resistance_evaluator, inertia_evaluator,
        resistance = std::vector<double>{},
        inertia = std::vector<double>{}](...) mutable
{
  const size_t element_count = airway_data.number_of_elements();
  resistance_evaluator(..., resize_scratch(resistance, element_count));
  inertia_evaluator(..., resize_scratch(inertia, element_count));
  ...
};
```

The scratch vectors are resized on each call, but capacity is reused for fixed-size model blocks. This avoids per-assembly vector allocation after the first call.

### Kelvin-Voigt Scratch Buffers

Kelvin-Voigt wall evaluators now reuse separate buffers for:

```text
resistance_derivative_q1
resistance_derivative_q2
inertia_derivative_q1
inertia_derivative_q2
viscous_wall_resistance_derivative_q1
viscous_wall_resistance_derivative_q2
```

The coefficient insertion order in `evaluate_tree_linearization_kelvin_voigt_wall(...)` is unchanged.

## Expected Counter Impact

The primary expected effect is lower allocation overhead in airway assembly callbacks.

Most relevant counter:

```text
tree_assembly_airways_s
```

Secondary effects may appear in full residual/Jacobian paths because those shared evaluator APIs were also changed.

This step is not expected to directly reduce:

```text
tree_assembly_solver_update_s
```

That cost still comes from solver coefficient-location resolution after `TreeLinearization` has been assembled.

## Verification

Built release targets:

```bash
cmake --build /scratch/Rodriguez/workspace/4C/4C/build/release \
  --target 4C benchmarktests_reduced_lung unittests_reduced_lung \
  --parallel 4
```

Ran reduced-lung unit tests:

```bash
ctest -R "^unittests_reduced_lung$" --output-on-failure
```

Result:

```text
100% tests passed, 0 tests failed out of 1
```

## Re-Measurement

From:

```text
/scratch/Rodriguez/workspace/4C/4C
```

run:

```bash
FOUR_C_REDUCED_LUNG_TREE_PROFILE=1 ./build/release/4C \
  ../files/reduced_lung_lung_tree_gen16_500steps_newton_tree.4C.yaml \
  ../output/gen16_newton_tree_profile
```

Compare `tree_assembly_airways_s` and `tree_assembly_s` against the Step 3 profile.
