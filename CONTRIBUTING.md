# Contributing to EnsembleQL

Contributions are welcome through GitHub issues and pull requests. Scientific behavior must remain explicit, testable, and documented.

## Development setup

Create an isolated Python environment, install an editable build, and run both test layers:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[test]"
python -m pytest -q
cmake -S . -B build -DENSEMBLEQL_BUILD_PYTHON=OFF -DENSEMBLEQL_BUILD_TESTS=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Use the commands in [`validation/README.md`](validation/README.md) when a change affects observable definitions, trajectory interpretation, units, periodic boundaries, or event extraction.

## Pull requests

- Keep changes focused and explain the scientific or user-facing motivation.
- Add tests for new behavior and boundary conditions.
- Update observable definitions and examples when semantics change.
- Do not commit generated validation trajectories.
- Keep public versions synchronized with `python tools/check_release.py`.
- Confirm the complete CI suite passes before requesting review.

By contributing, you agree that your contribution is licensed under the project's MIT License.
