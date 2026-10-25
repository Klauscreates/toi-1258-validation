"""Roman/Hubble microlensing architecture scaffold."""

from .candidates import Candidate, CandidateScore
from .config import ArchitectureConfig
from .mast_archive import query_jwst_proxima_observations, summarize_jwst_archive
from .proxima_mission import Mission2Config, build_mission_package, classify_habitability
from .simulation import paczynski_magnification, paczynski_flux

__all__ = [
    "ArchitectureConfig",
    "Candidate",
    "CandidateScore",
    "Mission2Config",
    "build_mission_package",
    "classify_habitability",
    "paczynski_flux",
    "paczynski_magnification",
    "query_jwst_proxima_observations",
    "summarize_jwst_archive",
]
