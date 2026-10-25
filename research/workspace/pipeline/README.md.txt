# Astro Candidate Pipeline

Saved: 2026-06-14

This is the reusable pipeline distilled from the Sirius / IRPCA-002 investigation.
Use it for future targets, planets, companion searches, archive hunts, or unusual
high-contrast candidates.

## Start Here

1. `FUTURE_THREAD_PROMPT.md`
2. `TARGET_MANIFEST_TEMPLATE.json`
3. `PLAYBOOK.md`
4. `PIPELINE_MAP.json`

## What This Pipeline Does

The pipeline turns a target hypothesis into a reproducible evidence package:

1. Define the target and predicted search geometry.
2. Convert separation/PA instructions into exact sky coordinates.
3. Query public archives.
4. Reduce high-contrast imaging.
5. Test sky-fixed versus detector-fixed behavior.
6. Run PCA/KLIP-style residual tests.
7. Run fake-source injection.
8. Check thermal, polarimetric, spectroscopic, and interferometric gates.
9. Record all negative evidence.
10. Package figures, tables, reports, scripts, and a paper draft.

## Multi-Messenger Survey Extension

The pipeline also includes a survey-scale branch for high-mass dynamical
anomalies. It combines observables without treating them as interchangeable:

- Gaia DR3 astrometric orbits rank massive unseen-companion proxies.
- SIMBAD identifies already classified systems.
- Chandra CSC 2.1 tests for a positional X-ray counterpart.
- TESS/MAST records time-domain coverage for later pixel-level analysis.
- Calibrated spectroscopy is a later chemistry gate, not an inference from
  broadband transit photometry.
- Euclid weak-lensing products belong to a separate extragalactic branch.

Run the pilot with:

```text
python scripts/multimessenger_governing_nodes_pilot.py --refresh-gaia
python scripts/tess_high_mass_proxy_pixel_parser.py
python scripts/hd85563_difference_image_localization.py
```

Current manifest:

`MULTIMESSENGER_GOVERNING_NODES_MANIFEST.json`

## Current Reference Implementation

The reference case is IRPCA-002:

`outputs/IRPCA002_CANDIDATE_PAPER_PACKAGE_2026-06-13/`

Key reference result:

```text
target: IRPCA-002
system: Sirius
separation: 1.185 arcsec
PA: 129.5 deg
east: +0.9146 arcsec
north: -0.7538 arcsec
SPHERE sky-minus max SNR: 16.8869
detector-fixed max SNR: 2.0665
NICMOS same-slot median SNR: 8.8034
STIS repeated projected-slot residual: 4041.876 Angstrom
```

## Future Thread Shortcut

In a new Codex thread, say:

```text
Use our saved Astro Candidate Pipeline from:
.\pipelines\astro_candidate_pipeline

Create a target manifest for [TARGET] and execute the pipeline.
```

## Important Claim Boundary

The pipeline is designed to preserve strong results without inflating them.
It separates:

- detection/candidate evidence;
- artifact rejection;
- archive access boundaries;
- physical interpretation;
- final classification.
