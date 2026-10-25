from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass(frozen=True)
class MastQueryPlan:
    """A serializable plan for a MAST observation query.

    This object deliberately stores query intent instead of performing network
    access at import time. Production workers can execute it with astroquery or
    the MAST API once credentials, rate limits, and storage policy are known.
    """

    mission: str
    target_name: str | None = None
    coordinates: str | None = None
    radius_degrees: float | None = None
    filters: tuple[str, ...] = ()
    proposal_ids: tuple[str, ...] = ()
    product_types: tuple[str, ...] = ("SCIENCE",)
    extra_criteria: dict[str, Any] = field(default_factory=dict)

    def to_astroquery_criteria(self) -> dict[str, Any]:
        criteria: dict[str, Any] = {"obs_collection": self.mission}
        if self.target_name:
            criteria["objectname"] = self.target_name
        if self.filters:
            criteria["filters"] = list(self.filters)
        if self.proposal_ids:
            criteria["proposal_id"] = list(self.proposal_ids)
        criteria.update(self.extra_criteria)
        return criteria


def m31_hubble_seed_plan() -> MastQueryPlan:
    """Return a conservative seed query for public HST M31 exploration."""

    return MastQueryPlan(
        mission="HST",
        target_name="M31",
        product_types=("SCIENCE",),
        extra_criteria={
            "dataproduct_type": "image",
            "intentType": "science",
        },
    )


def roman_future_wfi_plan(field_name: str) -> MastQueryPlan:
    """Return a future Roman WFI query intent.

    The exact MAST query fields may change before Roman operations, so this is
    a provenance-friendly placeholder rather than a hard-coded API contract.
    """

    return MastQueryPlan(
        mission="ROMAN",
        target_name=field_name,
        product_types=("SCIENCE",),
        extra_criteria={
            "instrument_name": "WFI",
            "dataproduct_type": "image",
            "calibration_level": [2, 3, 4],
        },
    )
