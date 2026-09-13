def calculate_gc_skew(sequence, window=1000, step=100):
    """
    Calculate GC skew across a DNA sequence using sliding windows.

    GC skew = (G - C) / (G + C)

    Args:
        sequence (str): DNA sequence.
        window (int): Size of each sliding window.
        step (int): Number of bases to move between windows.

    Returns:
        list[dict]: Window positions and GC-skew values.
    """
    sequence = sequence.upper()

    if window <= 0 or step <= 0:
        raise ValueError("window and step must be positive")

    if window > len(sequence):
        return []

    regions = []

    for start in range(0, len(sequence) - window + 1, step):
        end = start + window
        window_sequence = sequence[start:end]

        g_count = window_sequence.count("G")
        c_count = window_sequence.count("C")

        if g_count + c_count == 0:
            skew = 0.0
        else:
            skew = (g_count - c_count) / (g_count + c_count)

        regions.append({
            "start": start,
            "end": end,
            "gc_skew": skew
        })

    return regions


def calculate_cumulative_gc_skew(sequence):
    """
    Calculate cumulative GC skew across a DNA sequence.

    G increases the cumulative value by 1.
    C decreases the cumulative value by 1.

    Args:
        sequence (str): DNA sequence.

    Returns:
        list[float]: Cumulative GC-skew values.
    """
    sequence = sequence.upper()

    cumulative_skew = []
    current_skew = 0.0

    for base in sequence:
        if base == "G":
            current_skew += 1
        elif base == "C":
            current_skew -= 1

        cumulative_skew.append(current_skew)

    return cumulative_skew


def find_skew_regions(skew, cumulative_skew):
    """
    Identify candidate regions from GC-skew and cumulative-skew data.

    Args:
        skew (list[dict]): Sliding-window GC-skew results.
        cumulative_skew (list[float]): Cumulative GC-skew values.

    Returns:
        list[dict]: Candidate GC-skew regions.
    """
    if not skew:
        return []

    regions = []

    # Find minimum cumulative GC skew.
    min_position = cumulative_skew.index(min(cumulative_skew))

    for region in skew:
        start = region["start"]
        end = region["end"]

        # Regions near the cumulative-skew minimum
        # receive stronger consideration as OriC candidates.
        if start <= min_position <= end:
            regions.append({
                "start": start,
                "end": end,
                "gc_skew": region["gc_skew"],
                "cumulative_min_position": min_position
            })

    return regions
