# Implementation NextStep 9

## Goal

Create large `Reduced_Lung` runtime inputs from the external generation-10 lung tree file so the existing reduced-lung nonlinear solver workflows can be tested on the same geometry:

```text
/scratch/Rodriguez/workspace/4C/files/lung_tree_gen10.json
```

The source JSON cannot be used directly by 4C because it uses the external schema:

```text
nodepos
airways
alveolarclusters
```

The `Reduced_Lung` input parser expects topology and boundary-condition fields such as:

```text
node_coordinates
radius
bc_node_id
bc_function_id
```

## Added Field Data

Added:

```text
tests/input_files/reduced_lung_lung_tree_gen10_airways_fields.json
```

This file is an airways-only conversion of `lung_tree_gen10.json`.

Conversion summary:

- Original airway-tree node ids were remapped to dense 1-based 4C input node ids.
- Original airway entries were remapped to dense 1-based 4C input element ids.
- `node_coordinates` came from `nodepos`.
- `element_nodes`, `generation`, and `radius` came from `airways`.
- The root node is assigned one inlet pressure boundary condition using function 1.
- Every airway leaf outlet is assigned a pressure boundary condition using function 2.
- `alveolarclusters` were intentionally not converted in this first version.

Converted problem size:

```text
Nodes:                1024
Airway elements:      1023
Outlet leaf nodes:     512
Boundary conditions:   513
```

## Added Solver Inputs

Added three `Reduced_Lung` runtime YAML files that all reference the generated field JSON:

```text
tests/input_files/reduced_lung_lung_tree_gen10_nox.4C.yaml
tests/input_files/reduced_lung_lung_tree_gen10_newton_sparse.4C.yaml
tests/input_files/reduced_lung_lung_tree_gen10_newton_tree.4C.yaml
```

The files differ only in the nonlinear solver selection:

```yaml
nonlinear_solver: Nox
nonlinear_solver: NewtonSparse
nonlinear_solver: NewtonTree
```

All three use:

```yaml
SOLVER: "UMFPACK"
number_of_steps: 1
results_every: 1000000
```

The high `results_every` value suppresses routine VTK output for the one-step timing input, so manual runs focus on solver and assembly work rather than output writing.

## Test Registration

Registered the new YAML files in:

```text
tests/list_of_tests.cmake
```

Added CTest entries:

```cmake
four_c_test(TEST_FILE reduced_lung_lung_tree_gen10_nox.4C.yaml TIMEOUT 240)
four_c_test(TEST_FILE reduced_lung_lung_tree_gen10_newton_sparse.4C.yaml TIMEOUT 240)
four_c_test(TEST_FILE reduced_lung_lung_tree_gen10_newton_tree.4C.yaml TIMEOUT 240)
four_c_test(TEST_FILE reduced_lung_lung_tree_gen10_newton_tree.4C.yaml NP 2 TIMEOUT 240)
```

The serial entries compare `Nox`, `NewtonSparse`, and serial `NewtonTree`. The NP2 entry exercises the distributed `NewtonTree` runtime path with the same YAML file.

## Manual Run Commands

From `build/release`:

```text
./4C ../../tests/input_files/reduced_lung_lung_tree_gen10_nox.4C.yaml ../../../output_gen10_nox/
./4C ../../tests/input_files/reduced_lung_lung_tree_gen10_newton_sparse.4C.yaml ../../../output_gen10_sparse/
./4C ../../tests/input_files/reduced_lung_lung_tree_gen10_newton_tree.4C.yaml ../../../output_gen10_tree_serial/
mpirun -np 2 ./4C ../../tests/input_files/reduced_lung_lung_tree_gen10_newton_tree.4C.yaml ../../../output_gen10_tree_np2/
```

From the repository root:

```text
./build/release/4C tests/input_files/reduced_lung_lung_tree_gen10_nox.4C.yaml output_gen10_nox/
./build/release/4C tests/input_files/reduced_lung_lung_tree_gen10_newton_sparse.4C.yaml output_gen10_sparse/
./build/release/4C tests/input_files/reduced_lung_lung_tree_gen10_newton_tree.4C.yaml output_gen10_tree_serial/
mpirun -np 2 ./build/release/4C tests/input_files/reduced_lung_lung_tree_gen10_newton_tree.4C.yaml output_gen10_tree_np2/
```

## Current Limitations

- The conversion is airways-only; it does not convert `alveolarclusters` into terminal-unit elements.
- The generated case uses pressure boundary conditions at all terminal airway outlets.
- Terminal-unit volume from `alveolarclusters.volume` is not consumed by the current `Reduced_Lung` terminal-unit input path, so adding terminal units should be handled as a separate follow-up.
- The inputs are meant for solver-path timing and runtime stress, not for validating a physiological result against a reference.
