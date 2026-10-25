# Reproduction Guide

## Supported: bundled evidence audit

Requires Python 3.11 or newer. Run from this repository's root:

```shell
python -m venv .venv
```

Activate with `.venv\Scripts\Activate.ps1` on PowerShell, or
`source .venv/bin/activate` on macOS/Linux. Activation is optional if you
invoke the virtual environment's Python directly.

```shell
python -m pip install -e .
python -m toi1258
python -m unittest discover -s tests -v
python tools/validate_repository.py
```

The audit verifies the provenance checksums before calculation. It writes
`results/summary.json`, `results/audit.csv`, and `results/REPORT.md`.
There are no network calls, credentials, or runtime dependencies for this tier.
The evidence tables are real saved measurements, not simulated demonstration data.

To run from another directory, pass `--root /path/to/repository`.
Use `--output /path/to/results` to isolate output. The installed wheel contains
code, not the evidence archive; retain the repository checkout as the data root.

## Figures and development

```shell
python -m pip install -e ".[dev,figures]"
python tools/plot_overview.py
ruff check src tests tools
cffconvert --validate
python -m build
```

The figure is regenerated from the same bundled tables. It is not an illustration
of a confirmed physical object. CI is configured for Windows and Linux on Python
3.11 and 3.12; a configured matrix is not a claim that remote CI has already run.

## Not yet supported end-to-end: raw-data reconstruction

[Historical source](../research/README.md) records the original reduction logic.
It requires external archives, calibration products, and domain software not
bundled here. It is not executed by CI. Porting source paths alone does not make
a scientifically validated raw-data pipeline.

See [the data guide](DATA.md) for input products, units, and the work needed to
complete raw-data reproducibility. Never describe this offline audit as a fresh
DOLPHOT, SPHERE, or TRICERATOPS run.

