# Exact Input Reproduction

This adds an executable core reduction path to the offline evidence audit.
It is not a claim that every exploratory script has been ported or that any
candidate is confirmed. Read `reproduce/VERIFICATION.json` for executed coverage.

## Setup and retrieval

Use a checkout, not just the Python wheel. Commands below run at its root.
Python 3.13 on Windows was used for the core FITS reruns. The offline audit CI
is separate and does not download these large inputs. Install GitHub CLI and
authenticate to the private repository with your own account; no tokens are bundled.

```shell
python -m venv .venv-raw
```

Activate `.venv-raw/Scripts/Activate.ps1` on PowerShell or
`source .venv-raw/bin/activate` on a POSIX shell, then:

```shell
python -m pip install -r reproduce/requirements-raw.txt
python tools/raw_inputs.py --root ../external-data --release-assets --download --report ../input-verification.json
python tools/raw_inputs.py --root ../external-data
```

`reproduce/inputs.json` contains exact byte lengths, SHA-256 hashes, archive
identifiers, retrieval URLs and private release-asset mappings. Files are never
accepted on filenames alone. Changed upstream versions cause failure, not silent
substitution. Retained historical snapshots are provided where required. Inputs
stay outside Git; allow space for archives, extraction and intermediate products.
Interrupted extraction may need a new empty data directory; existing files are
not overwritten. `--select` and `--limit` support bounded retrieval checks.

The internal `--local-map` option was used to verify retained files without
duplicating the original cache. That machine-specific map is intentionally not
published. External users use the manifest-relative layout created by retrieval.

## Provenance and boundaries

Source and historical parameter choices are preserved under `research/`.
Wrappers substitute paths and executable locations, not evidence claims.
The manifest covers the inputs for the documented core paths, not an inventory
of every telescope archive ever consulted. Run logs write actual versions and
outputs into a new external directory. Each rerun refuses an existing output.

Reproduction of an estimator is not proof that its uncertainty model or physical
interpretation is correct. Historical manuscripts remain historical artifacts.

## TOI FITS and posterior paths

29 exact inputs total approximately 606 MB: 12 SPOC light curves, 12 target-pixel
products and five follow-up FITS. SPOC products are calibrated pipeline products,
not detector-raw telemetry. Follow-up snapshots are in the private data release.

```shell
python tools/run_raw.py --data-root ../external-data --output ../toi-fits-rerun
```

This stitches all twelve sectors, repeats the 6000-period BLS grid, phase-folded
metrics and target-pixel centroid checks. Outputs are in `evidence/tables/` under
the selected output. `reproduce/target.json` freezes the historical ephemeris.

Use a separate Python 3.12 environment for the historical TRICERATOPS sampler:

```shell
python -m pip install -r reproduce/requirements-triceratops.txt
python tools/run_triceratops.py --output ../toi-posterior-rerun
```

To explicitly connect the regenerated FITS light curve to posterior sampling:

```shell
python tools/run_triceratops.py --output ../toi-connected-posterior --binned-lightcurve ../toi-fits-rerun/evidence/tables/triceratops_folded_binned_lightcurve.csv
```

Defaults are five seeds (1258-1262), one million samples each. Compatibility
patches, the actual historical environment inventory, TRILEGAL realization,
contrast curves, binned light curve and star/flux-ratio/depth table are preserved.
The star table is a frozen input: this does not rerun current catalogs or the
network constructor's aperture contamination calculation. Sampling reproduction
is conditional on those historical inputs, not a fresh validation or measured mass.
The FITS wrapper also reruns the historical 25,000-row downsampling and folded
binning step. See the verification report for numerical comparison to the frozen
input. The five-seed verification uses the preserved historical table; the optional
argument selects a regenerated one and records its SHA-256 in the sampling log.
