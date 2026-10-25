# Scientific Status

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

- Reproduce the complete likelihood run with a reviewed software environment.
- Verify light-curve detrending, contamination, centroid and imaging assumptions.
- Audit the historical Python compatibility patch used for TRICERATOPS.
- Refresh the literature and archive status with dated queries.
- Seek independent mass constraints where available.

Methods: [TRICERATOPS](https://github.com/stevengiacalone/triceratops),
[TESS archive](https://archive.stsci.edu/missions-and-data/tess),
and [ExoFOP-TESS](https://exofop.ipac.caltech.edu/tess/).

