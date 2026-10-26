# Scientific Status

## Current result (corrected rerun, 2026-10-25)

The historical probabilities below are **superseded** by the
[corrected rerun](../rerun/README.md). The historical Gemini contrast curves
were corrupted: their separation column was overwritten with Δmag, so the run
effectively had no imaging constraint. The historical TRICERATOPS module also
predates the upstream 2026 likelihood and prior fixes.

The rerun uses current upstream TRICERATOPS, corrected contrast curves and a
freshly folded SPOC PDCSAP light curve. It gives FPP < 5×10⁻⁶ and NFPP ≈ 10⁻¹⁵
in every seed, both with public data only (TESS plus Gemini 2020) and with the
2026 Palomar AO curve added. The SPOC multi-sector centroid excludes
TIC 288132263 at 5.8σ, and per-sector difference images exclude it at 8.5σ.
PDCSAP per-sector depths are consistent (792 ± 13 ppm). These are statistical
validation results, not an RV mass, and they are conditional on the scenario
set, the priors and the inputs described in the rerun README.

## Historical evidence

The retained analysis snapshot reports a twelve-sector TESS signal near
10.261376896 days and a BLS depth of 722.467 ppm. Five saved million-sample
TRICERATOPS seeds report FPP approximately 2.37e-4 to 6.92e-4 and NFPP
approximately 4.26e-5 to 6.50e-5.

The offline audit reconstructs those probabilities from the saved scenario
tables. It does not rerun the likelihood model, verify all input priors, establish
novelty, or check today's archive classification. Low probabilities are
conditional on the scenarios, data preparation, contrast curves, and priors.

## Validation is not dynamical confirmation

No RV mass or atmospheric spectrum is measured by this repository. A substantial
RV campaign listed in metadata is not itself a public reduced RV time series.
No assumption is made that the July 2026 catalog classification is current.

The percentage difference between two period estimates is not a calibrated
period uncertainty. Different transit-depth estimators are kept separate.
Centroid non-detection of an offset constrains blends; it does not by itself
prove that every false-positive scenario has been excluded.

## Reproduced here

FPP is the sum of scenarios other than TP, PTP and DTP; NFPP is the sum of
NTP, NEB and NEBx2P scenarios. Probabilities must be finite, in [0, 1], and sum
to unity. The recomputed values must agree with the saved full-scenario sweep.
The 0.001 thresholds reproduce the historical report and are not a universal
guarantee of journal acceptance.

## Before an independent validation claim

- ~~Reproduce the complete likelihood run with a reviewed software environment.~~
  Done in the [corrected rerun](../rerun/README.md), with upstream TRICERATOPS
  pinned to a commit.
- ~~Verify light-curve detrending, contamination, centroid and imaging assumptions.~~
  Done there. This check found the contrast-curve corruption.
- ~~Audit the historical Python compatibility patch used for TRICERATOPS.~~
  It was an older upstream module, replaced by current upstream.
- Refresh the literature and archive status with dated queries. As of
  2026-10-25 it is TFOPWG "PC", and no publication was found.
- Coordinate with the TFOP data owners before publishing (ExoFOP data-use policy).
- Seek independent mass constraints where available.

Methods: [TRICERATOPS](https://github.com/stevengiacalone/triceratops),
[TESS archive](https://archive.stsci.edu/missions-and-data/tess),
and [ExoFOP-TESS](https://exofop.ipac.caltech.edu/tess/).

