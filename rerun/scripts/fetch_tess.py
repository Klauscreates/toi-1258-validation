"""Download SPOC 120-s light curves and target pixel files for TOI-1258."""
import sys
from pathlib import Path
import lightkurve as lk

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
for kind, search in (("lc", lk.search_lightcurve), ("tpf", lk.search_targetpixelfile)):
    res = search("TIC 288132261", mission="TESS", author="SPOC", exptime=120)
    print(kind, len(res), "products", flush=True)
    res.download_all(download_dir=str(out / kind))
