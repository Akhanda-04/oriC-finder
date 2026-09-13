"""
Scoring functions for Ori-C Finder.

This module combines multiple features of candidate regions into
a single score and ranks candidate regions from highest to lowest.

Expected candidate dictionary structure:

{
    "start": 3920000,
    "end": 3925000,
    "length": 5000,
    "gc_skew": 0.18,
    "motif_count": 6
}

The scoring system is intentionally simple and interpretable.
"""


def score_candidate(
    candidate,
    motif_weight=2.0,
    skew_weight=3.0,
    length_weight=1.0
):
    """
    Calculate an Ori-C likelihood score for one candidate region.

    Parameters
    ----------
    candidate : dict
        Candidate region containing:
        - motif_count
        - gc_skew
        - length

    motif_weight : float
        Weight given to DnaA motif evidence.

    skew_weight : float
        Weight given to GC-skew evidence.

    length_weight : float
        Weight given to candidate length.

    Returns
    -------
    float
        Candidate score.
    """

    motif_count = candidate.get("motif_count", 0)
    gc_skew = abs(candidate.get("gc_skew", 0.0))
    length = candidate.get("length", 0)

    # Normalize length so very large regions do not dominate the score.
    length_score = min(length / 10000, 1.0)

    score = (
        motif_weight * motif_count
        + skew_weight * gc_skew
        + length_weight * length_score
    )

    return score


def rank_candidates(
    candidates,
    motif_weight=2.0,
    skew_weight=3.0,
    length_weight=1.0
):
    """
    Score and rank candidate Ori-C regions.

    Parameters
    ----------
    candidates : list of dict
        Candidate regions.

    Returns
    -------
    list of dict
        Candidates sorted from highest score to lowest score.
        Each candidate receives a new "score" field.
    """

    scored_candidates = []

    for candidate in candidates:
        candidate_copy = candidate.copy()

        candidate_copy["score"] = score_candidate(
            candidate_copy,
            motif_weight=motif_weight,
            skew_weight=skew_weight,
            length_weight=length_weight
        )

        scored_candidates.append(candidate_copy)

    # Highest-scoring candidate first
    scored_candidates.sort(
        key=lambda candidate: candidate["score"],
        reverse=True
    )

    return scored_candidates
