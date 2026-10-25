# Contributing

Start with [the evidence status](docs/SCIENCE.md) and [reproduction guide](docs/REPRODUCING.md).

1. Create a branch for a focused change.
2. Install the project with `python -m pip install -e ".[dev,figures]"`.
3. Run `python -m unittest discover -s tests -v`, `ruff check src tests tools`,
   and `python tools/validate_repository.py`.
4. Include commands, environment versions, units, controls, and expected outputs
   in the pull request. Do not submit large raw products or personal information.

Do not overwrite historical evidence to make a test pass. A changed measurement
needs a new provenance record, a justified method change, and updated scientific
interpretation. Test failures should not be resolved by removing unfavorable
measurements. Keep candidate, validation, and confirmation language distinct.

Until licensing is resolved, discuss proposed contributions before submitting
third-party material or code you do not own.

