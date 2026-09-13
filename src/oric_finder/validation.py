def calculate_overlap(predicted_start, predicted_end,
                      known_start, known_end):
    """
    Calculate the number of bases shared by predicted and known oriC.
    """

    overlap_start = max(predicted_start, known_start)
    overlap_end = min(predicted_end, known_end)

    if overlap_start >= overlap_end:
        return 0

    return overlap_end - overlap_start


def calculate_coverage(predicted_start, predicted_end,
                       known_start, known_end):
    """
    Calculate how much of the known oriC is covered
    by the predicted region.
    """

    overlap = calculate_overlap(
        predicted_start, predicted_end,
        known_start, known_end
    )

    known_length = known_end - known_start

    if known_length == 0:
        return 0.0

    return overlap / known_length


def calculate_distance(predicted_start, predicted_end,
                       known_start, known_end):
    """
    Calculate the distance between predicted and known oriC regions.

    Returns 0 if the regions overlap.
    """

    if predicted_end >= known_start and known_end >= predicted_start:
        return 0

    if predicted_end < known_start:
        return known_start - predicted_end

    return predicted_start - known_end


def validate_prediction(predicted_start, predicted_end,
                        known_start, known_end,
                        coverage_threshold=0.5,
                        distance_threshold=5000):
    """
    Determine whether an oriC prediction agrees with a known oriC.
    """

    overlap = calculate_overlap(
        predicted_start, predicted_end,
        known_start, known_end
    )

    coverage = calculate_coverage(
        predicted_start, predicted_end,
        known_start, known_end
    )

    distance = calculate_distance(
        predicted_start, predicted_end,
        known_start, known_end
    )

    validated = (
        coverage >= coverage_threshold
        or distance <= distance_threshold
    )

    return {
        "overlap": overlap,
        "coverage": coverage,
        "distance": distance,
        "validated": validated
    }
