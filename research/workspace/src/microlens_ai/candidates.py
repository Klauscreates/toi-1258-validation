from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class CandidateScore:
    """Ranking metrics produced by the model and physics-vetting layer."""

    microlensing_probability: float
    artifact_probability: float
    variable_probability: float
    uncertainty: float

    def vetting_priority(self) -> float:
        """Higher means more worth human review."""

        signal = self.microlensing_probability
        contamination = 0.5 * self.variable_probability + self.artifact_probability
        confidence_penalty = 0.25 * self.uncertainty
        return signal - contamination - confidence_penalty


@dataclass(frozen=True)
class Candidate:
    """Minimal candidate packet metadata.

    Production packets should also include archive product IDs, cutout URIs,
    fit summaries, and crossmatch tables.
    """

    candidate_id: str
    ra_degrees: float
    dec_degrees: float
    field_id: str
    first_mjd: float
    last_mjd: float
    bands: tuple[str, ...]
    score: CandidateScore | None = None
    provenance: dict[str, str] = field(default_factory=dict)

    @property
    def duration_days(self) -> float:
        return self.last_mjd - self.first_mjd
