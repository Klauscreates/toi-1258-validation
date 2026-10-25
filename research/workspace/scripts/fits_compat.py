from __future__ import annotations


def disable_broken_lzw_detection() -> None:
    """Avoid Astropy FITS crashes when a partial uncompresspy package is present."""
    try:
        from astropy.io.fits import file as fits_file

        if not getattr(fits_file, "HAS_UNCOMPRESSPY", False):
            return

        import uncompresspy

        if not hasattr(uncompresspy, "LZWFile"):
            fits_file.HAS_UNCOMPRESSPY = False
    except Exception:
        return
