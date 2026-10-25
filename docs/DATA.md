# Data and Provenance

This package curates the July 2026 publication evidence tables and June 2026
saved posterior scenario tables. The export manifest records exact source and
public-export SHA-256 hashes. Absolute machine paths are removed or reduced
to filenames without rounding scientific values.

| Files | Units or role |
| --- | --- |
| stitched_bls_summary.csv | Period/duration in days, depth in ppm; historical BLS estimate |
| stitched_transit_metrics.csv | Phase-folded primary/secondary and other diagnostic measurements |
| stitched_transit_depth_by_sector.csv | Sector-specific depth diagnostics |
| triceratops_folded_binned_lightcurve.csv | Time from mid-transit in days, normalized flux/error |
| triceratops_seed_sweep_N1000000_summary.csv | Seed identifiers, sample counts, historical probabilities |
| probabilities_full_seed_*.csv | Full scenario posteriors, not simulated replacement data |
| triceratops_contrast_curve_*.csv | Historical contrast curves; see column names for separation and delta magnitude |
| spoc_dv_centroid_summary.csv | SPOC offsets and errors in arcsec; multiple pipeline products |
| mast_tess_product_inventory.csv | Archived product identifiers for reacquisition |

## Raw reconstruction

Acquire TESS products using the inventory through [MAST](https://mast.stsci.edu/).
Contrast/follow-up inputs originate in [ExoFOP-TESS](https://exofop.ipac.caltech.edu/tess/).
The [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/) provides
catalog context, but no live status is asserted by this snapshot.

A full rerun requires original light curves and pixel products, detrending
settings, the complete TRICERATOPS environment and reviewed compatibility changes,
and the original imaging inputs. The historical modified likelihood source is
retained with its upstream MIT attribution under research/third_party/. It needs
the rest of its upstream package and a reviewed environment. Thus the offline
posterior check is reproducible now; the full historical posterior generation
is not yet packaged as an end-to-end supported command.

The portable historical SPOC XML parser may be inspected in research/. Its
output still needs instrument/scientific validation; passing a software test is
not independent confirmation of a planet.

