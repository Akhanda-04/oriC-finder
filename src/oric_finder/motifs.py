"""Functions for detecting and scoring DnaA-box motifs."""

from typing import List, Tuple


def score_motif(sequence: str, motif: str) -> int:
    """
    Score how well a DNA sequence matches a motif.

    The score is the number of matching positions.

    Parameters
    ----------
    sequence : str
        DNA sequence to evaluate.
    motif : str
        Reference motif.

    Returns
    -------
    int
        Number of matching positions.

    Raises
    ------
    ValueError
        If the sequences have different lengths.
    """
    sequence = sequence.upper()
    motif = motif.upper()

    if len(sequence) != len(motif):
        raise ValueError("Sequence and motif must have the same length.")

    return sum(base == reference for base, reference in zip(sequence, motif))


def find_dnaa_motifs(
    sequence: str,
    consensus: str = "TTATCCACA",
    min_score: int = 7,
) -> List[Tuple[int, str, int]]:
    """
    Find potential DnaA-box motifs in a DNA sequence.

    A sliding window is used to compare every region of the sequence
    against the DnaA-box consensus sequence.

    Parameters
    ----------
    sequence : str
        DNA sequence to scan.
    consensus : str, optional
        DnaA-box consensus motif. Default is TTATCCACA.
    min_score : int, optional
        Minimum number of matching positions required to report a motif.
        Default is 7.

    Returns
    -------
    list of tuple
        Each tuple contains:

        (position, motif, score)

        Position is zero-based.

    Raises
    ------
    ValueError
        If the sequence contains invalid DNA bases.
        If min_score is outside the valid range.
    """
    sequence = sequence.upper()
    consensus = consensus.upper()

    valid_bases = {"A", "T", "G", "C"}

    if any(base not in valid_bases for base in sequence):
        raise ValueError("Sequence contains invalid DNA bases.")

    if any(base not in valid_bases for base in consensus):
        raise ValueError("Consensus contains invalid DNA bases.")

    motif_length = len(consensus)

    if min_score < 0 or min_score > motif_length:
        raise ValueError(
            f"min_score must be between 0 and {motif_length}."
        )

    results = []

    for position in range(len(sequence) - motif_length + 1):
        window = sequence[position:position + motif_length]

        score = score_motif(window, consensus)

        if score >= min_score:
            results.append((position, window, score))

    return results
