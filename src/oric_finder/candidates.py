def identify_candidate_regions(
    skew_regions,
    motifs,
    sequence_length,
    min_length=500,
    max_length=5000
):
    """
    Identify candidate Ori-C regions using GC-skew regions
    and DnaA-box-like motifs.

    Parameters
    ----------
    skew_regions : list of dict
        Regions returned by find_skew_regions().

    motifs : list of tuple
        Motifs returned by find_dnaa_motifs().
        Each tuple is:
        (position, motif, score)

    sequence_length : int
        Length of the genome.

    min_length : int
        Minimum candidate region length.

    max_length : int
        Maximum candidate region length.

    Returns
    -------
    list of dict
        Candidate regions with GC-skew and motif evidence.
    """

    candidates = []

    for region in skew_regions:
        start = region["start"]
        end = region["end"]

        # Handle circular genome wrap-around
        if start <= end:
            length = end - start
        else:
            length = (sequence_length - start) + end

        if not (min_length <= length <= max_length):
            continue

        # Count DnaA motifs inside the candidate region
        motif_count = 0

        for position, motif, score in motifs:

            if start <= end:
                inside = start <= position < end
            else:
                inside = (
                    position >= start
                    or position < end
                )

            if inside:
                motif_count += 1

        candidates.append({
            "start": start,
            "end": end,
            "length": length,
            "gc_skew": region["gc_skew"],
            "motif_count": motif_count,
            "cumulative_min_position":
                region["cumulative_min_position"]
        })

    return candidates


def merge_regions(regions, max_gap=500):
    """
    Merge overlapping or nearby genomic regions.
    """

    if not regions:
        return []

    regions = sorted(
        regions,
        key=lambda region: region["start"]
    )

    merged = [regions[0].copy()]

    for region in regions[1:]:
        previous = merged[-1]

        if region["start"] <= previous["end"] + max_gap:

            previous["end"] = max(
                previous["end"],
                region["end"]
            )

            previous["length"] = (
                previous["end"] - previous["start"]
            )

            previous["motif_count"] += region["motif_count"]

            # Keep the strongest GC-skew evidence
            if abs(region["gc_skew"]) > abs(previous["gc_skew"]):
                previous["gc_skew"] = region["gc_skew"]

        else:
            merged.append(region.copy())

    return merged