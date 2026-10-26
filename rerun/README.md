# Corrected Rerun (2026-10-25)

This directory supersedes the historical TRICERATOPS probabilities in `data/`.
The historical files are kept unchanged for provenance. Their problems are
documented below and in [scientific status](../docs/SCIENCE.md).

## Why a rerun was needed

1. **Corrupted contrast curves.** In both historical Gemini/'Alopeke curves
   (`data/triceratops_contrast_curve_{562,832}nm.csv`, and the copies in
   `reproduce/historical-inputs/`), the separation column was overwritten with the
   Δmag column, so every row reads (Δm, Δm). The 0.5″ contrast used was 0.5 mag,
   while the measured values are 4.04 mag (562 nm) and 6.39 mag (832 nm). Only
   the 832 nm curve was passed to TRICERATOPS. In effect, the historical run had
   no high-resolution imaging constraint.
2. **Outdated TRICERATOPS.** The historical likelihood module predates the upstream
   2026 fixes. Those fixes changed the evidence calculation (the `exp(lnL + 600)`
   underflow), the background-companion prior log base, the BEB/NEB collision
   masks, and added support for multiple contrast curves.
3. **Depth estimate.** The historical stitched per-sector depths (517–694 ppm)
   are inconsistent with each other and biased low. SPOC PDCSAP depths measured
   with a local per-transit baseline are consistent with each other.

## Results

All numbers below are in [outputs/SUMMARY.json](outputs/SUMMARY.json), which
`scripts/summarize.py` regenerates from the committed result files.

**Transit depth (12 sectors, S14–S74):** PDCSAP 792 ± 13 ppm, per-sector
χ² = 16.4 for 11 degrees of freedom.

**Transit fit** ([outputs/transit_fit.json](outputs/transit_fit.json), `scripts/fit_transit.py`):
this is an emcee fit of the 25 transits with two-sided baselines (17,097 points). The free
parameters are T0, P, Rp/R★, b, ρ★, the Kipping q1 and q2, and jitter. The orbit is circular.
The run used 40 walkers for 60,000 steps; the longest autocorrelation time is 651 steps, so
the post-burn chain spans at least 45 autocorrelation times. Thinned posterior samples are
in `outputs/transit_fit.chain.npy`.

- P = 10.261314 ± 0.000011 d, and T0 = 2459430.6513 ± 0.0005 BJD_TDB.
- Rp/R★ = 0.0273 (+0.0005, −0.0004) and b = 0.32 (+0.23, −0.22).
- Rp = 4.44 ± 0.22 R⊕, using R★ = 1.49 ± 0.07 R☉ from TIC v8.
- The fitted ρ★ = 0.56 (+0.09, −0.18) g cm⁻³ agrees with the TIC density of the target
  (0.69 ± 0.15) at −0.7σ, which is evidence that the transit is on the target.

**Odd/even and secondary** ([outputs/vetting_checks.json](outputs/vetting_checks.json)):
the odd transits are 784 ± 17 ppm and the even transits 796 ± 18 ppm, a 0.5σ difference.
At phase 0.5 the flux change is −21 ± 12 ppm, so any secondary eclipse is below 35 ppm (3σ).

**Figure** for the Research Note: `scripts/make_figure.py`.

**Depth versus aperture:** from a 2 px to a 4 px radius, the depth gets
shallower in 6 of 10 usable sectors (median change −24 ppm). That is the
dilution pattern expected for an on-target source. Only S74 gets deeper
(+53 ppm, about 1σ). S60 and S73 have only one or two transits each and
strong SAP systematics, so they were excluded. ExoFOP's "possible
depth-aperture correlation" note comes from a QLP S74 reduction.

**Centroids and TIC 288132263** (T2: 17.0″, PA 64°, ΔT = 4.4; it would need
a 4.6% eclipse to be the source):

- The SPOC DV multi-sector solution (S14–S78) is 0.2σ from the target and
  5.8σ from T2.
- This rerun's per-sector difference images (10 sectors with at least two
  transits) give a mean offset of (+3.0″, −1.3″) ± 1.8″ per axis, which
  excludes T2 at 8.5σ.
- Single-sector SPOC outliers: S26 points south, away from T2, and still
  excludes T2 at 8.0σ. S74, the only camera-2 sector, is 4.9σ from the target
  and 3.5σ from T2, so it is inconsistent with both and treated as a
  sector-specific systematic.

**TRICERATOPS** (upstream commit `d20875c27d`, N = 10⁶ per seed, folded
PDCSAP light curve in 10-minute bins):

| Configuration | Seeds | max FPP | max NFPP |
|---|---|---|---|
| Gemini 562 + 832 nm, Gaia DR3 background (public data only) | 5 | 4.7×10⁻⁶ | 1.1×10⁻¹⁵ |
| Gemini + Palomar Kcont, Gaia DR3 background | 5 | 4.7×10⁻⁶ | 1.1×10⁻¹⁵ |
| Gemini 562 + 832 nm, frozen historical TRILEGAL | 3 | 1.6×10⁻¹⁰ | 5.1×10⁻¹⁵ |
| Gemini + Palomar Kcont, frozen historical TRILEGAL | 3 | 1.6×10⁻¹⁰ | 5.1×10⁻¹⁵ |
| No contrast curves (sensitivity check) | 3 | 2.6×10⁻⁴ | 9.9×10⁻¹⁶ |

- TP carries 94–98% of the probability, PTP 2–6%, and DTP below 0.1%. All
  three put the planet on the target star.
- When the FPP is this small it varies a lot between seeds, so quote it as an
  upper bound (FPP < 10⁻⁵), not a point value.
- The no-curve run reproduces the historical FPP range, which is consistent
  with the historical curves having provided no constraint.

**LCO 2022-05-08 (zs).** The observation starts after ingress, so it covers
the second half of the transit and egress. A simple trapezoid fit gives
1.5 ± 0.7 ppt on target, a marginal detection, so it is not used as
independent on-target confirmation. A box fit to T2's light curve shows an
apparent 21 ppt "dip", but T2 brightens smoothly through the night while the
sky background halves ([outputs/lco_T1_T2.png](outputs/lco_T1_T2.png)). That
is sky-subtraction error on a faint star in egress-only coverage, not an
eclipse. KeplerCam (2021) saturated the target. Upstream TRICERATOPS has no
external light-curve input. The TRICERATOPS+ fork does, but it predates the
2026 fixes, so it was not used.

## Limits

- Upstream lists the Gaia DR3 background population as still awaiting
  validation. The frozen-TRILEGAL runs agree with it here.
- The contrast-curve filters were mapped as 562 nm → Vis, 832 nm → TESS,
  Kcont → K.
- No RV mass is measured. Statistical validation is not dynamical confirmation.
- The Palomar AO curve (2026-04-07, PI Ciardi) may still be inside the
  12-month TFOP proprietary period. It is **not redistributed here**.
  `prep_inputs.py` reads it from an ExoFOP download. The result does not
  depend on it.
- The Gemini curves (2020-02-15, PI Crossfield, contact S. Howell) are
  included in corrected form. Cite them through ExoFOP as its data-use
  policy describes.

## Reproduce

Use a separate environment ([requirements.txt](requirements.txt)) and keep
downloads outside the repository:

```shell
python rerun/scripts/fetch_tess.py ../toi1258-data/tess
python rerun/scripts/prep_inputs.py ../toi1258-data data ../toi1258-out/inputs
python rerun/scripts/depth_aperture.py ../toi1258-data/tess ../toi1258-out/depth_aperture.json
python rerun/scripts/diff_centroids.py ../toi1258-data/tess ../toi1258-out/diff_centroids.json
python rerun/scripts/lco_check.py ../toi1258-data/exofop/TIC288132261-01_20220508_LCO-McD-1m0_zs_measurements.tbl ../toi1258-out/lco_check.json
python rerun/scripts/run_fpp.py --data ../toi1258-data --inputs ../toi1258-out/inputs --out ../toi1258-out/fpp_gemini --curves gemini
python rerun/scripts/fit_transit.py ../toi1258-data rerun/outputs/transit_fit.json --steps 60000
python rerun/scripts/vetting_checks.py ../toi1258-data rerun/outputs/transit_fit.json rerun/outputs/vetting_checks.json
python rerun/scripts/summarize.py rerun/outputs data/spoc_dv_centroid_summary.csv
python rerun/scripts/make_figure.py ../toi1258-data . ../toi1258-out/figure.png
```

`prep_inputs.py` expects these ExoFOP files in `../toi1258-data/exofop/`:
`TOI1258I-ef20200215-{562,832}_sensitivity.dat` and
`TOI1258I-dc20260407-Kcont_plot.tbl`. The LCO table is needed only for
`lco_check.py`. All four are on the [target's ExoFOP page](https://exofop.ipac.caltech.edu/tess/target.php?id=288132261)
(file IDs 237544, 237543, 1108115 and 880127), downloadable as
`https://exofop.ipac.caltech.edu/tess/get_file.php?id=<ID>`. The `--curves` options are `gemini`, `all` (adds Palomar),
and `none`. Add `--trilegal reproduce/historical-inputs/288132261_TRILEGAL.csv`
for the frozen-population comparison.

The TESS light curves and pixel files (about 560 MB) come from MAST through
lightkurve. TRICERATOPS queries the TIC and Gaia at setup.
