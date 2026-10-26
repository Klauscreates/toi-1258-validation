# Changelog

## 0.3.0 - 2026-10-25

- Corrected rerun under `rerun/`, with scripts, corrected Gemini contrast
  curves, a folded PDCSAP light curve, and all result tables plus a
  regenerable `SUMMARY.json`.
- Found that both historical contrast-curve CSVs had their separation column
  overwritten with Δmag. They are documented and kept unchanged for provenance.
- Reran TRICERATOPS at upstream commit d20875c27d: FPP < 5×10⁻⁶, NFPP ≈ 10⁻¹⁵,
  with public data only and with the 2026 Palomar curve, using Gaia DR3 and
  frozen TRILEGAL backgrounds.
- Added depth-versus-aperture, difference-image centroid (TIC 288132263
  excluded at 8.5σ) and LCO light-curve checks.
- Historical probabilities are superseded. The offline audit of the historical
  tables is unchanged.

## 0.2.0 - 2026-10-25

- Exact input manifests and protected data-release retrieval.
- Historical FITS reduction wrappers, tested environments and coverage reports.
- PHAT source recovery and frozen-input TOI posterior resampling where applicable.


## 0.1.2 - 2026-10-25

- Compare regenerated floating-point summaries with a strict 1e-12 relative
  tolerance across platforms; provenance hashes remain exact.
- No changes to inputs, probability calculations, or scientific conclusions.

## 0.1.1 - 2026-10-25

- Expanded project scripts, shared library, historical versions and manuscript sources.
- Added source coverage checks and explicit records for empty original files.
- Preserved LF text across operating systems so provenance hashes survive checkout.
- Published private GitHub repositories; no scientific confirmation claim is implied.

## 0.1.0 - 2026-10-25

- Curated historical measurements and reference source into a standalone project.
- Added offline evidence audit, input validation, checksums, tests, and CI.
- Documented scientific limitations separately from numerical results.
- Added a reproducible overview figure and citation metadata.
- No new archive reduction or publication acceptance is claimed by this release.

