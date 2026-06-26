# 4C Opencode Instructions

## Project Context

4C is a C++20 CMake/Ninja multiphysics research code with Python utilities and Sphinx/Doxygen documentation. Most production code lives under `src`, application entry points under `apps`, end-to-end input tests under `tests`, legacy unit tests under `unittests`, and developer tooling under `utilities`.

## Working Rules

- Keep changes small and aligned with nearby module patterns.
- Read the relevant source, CMake, and tests before editing; this codebase has legacy and newer module layouts side by side.
- Do not edit generated or local build artifacts under `build`, `.cache`, `utilities/python-venv`, or other ignored output directories.
- Do not edit `CMakeUserPresets.json` unless explicitly asked; it is local, ignored, and machine-specific.
- Add or update tests for behavior changes when practical. If verification is not practical because dependencies or a build are unavailable, state that clearly.
- Prefer ASCII in code and scripts. The pre-commit setup checks non-ASCII characters for C/C++, Python, CMake, and shell files.

## C++ Conventions

- Use C++20 and follow the C++ Core Guidelines where they fit the existing code.
- Add Doxygen comments for public C++ interfaces when introducing or changing API behavior.
- Naming conventions: namespaces, classes, structs, and enums use `CamelCase`; functions use `snake_case`; variables use `snake_case` or lower camel case; private members end with `_`; define flags use `FOUR_C_*`.
- Avoid single-letter variable names except conventional loop indices.
- Prefer `enum class` for new enums.
- Avoid new define flags. Use `FOUR_C_ENABLE_ASSERTIONS` for debug-only checks when needed.
- Avoid header-in-header inclusion where practical; prefer forward declarations, especially via `.fwd.hpp` files for repeated or external-library declarations.
- Pass parameters by `const` reference by default. Use `std::unique_ptr` or `std::shared_ptr` only when ownership or lifetime management requires it.
- Prefer specific parameter structs or container classes over `Teuchos::ParameterList` for new code.
- Keep code const-correct, especially function parameters and member functions.

## CMake Conventions

- Use out-of-source builds only. Never configure CMake in the source directory.
- Use presets when configuring. From the source root, list presets with `cmake . --list-presets`.
- This checkout currently has local presets `debug` and `release` in `CMakeUserPresets.json`; they build under `build/debug` and `build/release`.
- Configure from the source root with `cmake --preset=debug` or `cmake --preset=release` when those local presets exist.
- Build with `cmake --build build/debug --parallel <jobs>` or `ninja -C build/debug -j <jobs>` after configuring the `debug` preset.
- In CMake code, prefer target-based modern CMake. Project variables start with `FOUR_C_`; avoid modifying `CMAKE_` variables inside project CMake files unless there is a clear existing pattern.
- Project CMake helper functions start with `four_c_` in code and are documented as the `FOUR_C_` family in developer docs.

## Development Setup

- For development hooks and Python tooling, run `./utilities/set_up_dev_env.sh [optional-python]` once from the source root. It requires Python 3.12 or newer, creates `utilities/python-venv`, installs the editable Python package with development extras, and installs pre-commit hooks.
- After setup, run targeted hooks with `utilities/python-venv/bin/pre-commit run --files <changed files>`.
- Run all hooks with `utilities/python-venv/bin/pre-commit run --all-files` only when warranted; it can be expensive.
- Formatting is managed by pre-commit: clang-format for C/C++, cmake-format for CMake, black for Python, yamlfmt for YAML, and project-specific checks for filenames, includes, header guards, input files, and commit messages.

## Testing

- Run tests from the configured build directory with `ctest`.
- Use `ctest -R <regex>` for targeted tests, `ctest -R unittests` for unit tests, and `ctest -L minimal` for the minimal suite.
- For a single input file under `tests/input_files`, use `ctest -R <input_file>` with the file name or a specific regex.
- Unit tests usually live next to module source under `src/<module>/tests` and use `*_test.cpp` names. They are picked up by CMake through `four_c_auto_define_tests()`.
- Legacy unit tests live under `unittests`. For new legacy-style unit tests, mirror the tested file path, add the test to the local `CMakeLists.txt`, include GoogleTest, open an anonymous namespace, and follow nearby examples.
- In unit test executables, `FOUR_C_THROW` is replaced by a throwing version that raises `Core::Exception`; test errors with `EXPECT_THROW(<code>, Core::Exception)`.
- If adding input files under `tests/input_files`, make sure they are represented in `tests/list_of_tests.cmake` unless an existing exception applies.

## Documentation

- General documentation sources are under `doc/documentation`; source API docs use Doxygen.
- Documentation builds require the relevant CMake options. Build targets include `documentation` and `doxygen`, for example `cmake --build <build-dir> --target documentation`.
