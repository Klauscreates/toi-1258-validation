from __future__ import annotations

from .candidates import Candidate, CandidateScore


def rank_candidates(candidates: list[Candidate]) -> list[Candidate]:
    """Sort candidates by human-vetting priority."""

    def priority(candidate: Candidate) -> float:
        if candidate.score is None:
            return float("-inf")
        return candidate.score.vetting_priority()

    return sorted(candidates, key=priority, reverse=True)


def calibrated_score(
    microlensing_logit: float,
    variable_logit: float,
    artifact_logit: float,
    uncertainty: float,
) -> CandidateScore:
    """Create a simple normalized score from three class logits."""

    probs = _softmax([microlensing_logit, variable_logit, artifact_logit])
    return CandidateScore(
        microlensing_probability=probs[0],
        variable_probability=probs[1],
        artifact_probability=probs[2],
        uncertainty=max(0.0, min(1.0, uncertainty)),
    )


def _softmax(values: list[float]) -> list[float]:
    shifted = [value - max(values) for value in values]
    exps = [pow(2.718281828459045, value) for value in shifted]
    total = sum(exps)
    return [value / total for value in exps]
