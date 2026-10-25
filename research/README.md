# Historical Analysis Source

## Expanded Source Collection

`workspace/scripts/` contains project scripts and their local import dependencies.
`workspace/src/` includes the shared microlens_ai package. Original tests,
configuration, protocols and notes are retained alongside it. `snapshots/`
preserves distinct archived source versions, deduplicated by original hash.

[SOURCE_INVENTORY.json](SOURCE_INVENTORY.json) records origins, hashes and exclusions.
Empty source files remain empty and are flagged, not filled with invented code.
Historical conclusions may be superseded; current scientific status takes precedence.

These source snapshots preserve the original analysis methods for inspection.
They are **not supported end-to-end entry points**. Machine-specific prefixes
were replaced with relative external-data locations. Original and exported
SHA-256 hashes and transformations are recorded in [the manifest](../data/manifest.json).

Historical comments or report strings may contain interpretations that have
since changed. [Current scientific status](../docs/SCIENCE.md) takes precedence.
No external executables, calibration libraries, or third-party patched source
are redistributed except the historical modified TRICERATOPS likelihood module,
included with its upstream MIT license. Do not run these scripts against the
working evidence folder.

The tested public entry point is documented in [reproduction](../docs/REPRODUCING.md).
Making raw reduction portable requires explicit input manifests, calibrated
software environments, runtime tests, and review of output units and fit quality.

