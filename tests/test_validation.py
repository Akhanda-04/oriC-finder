from oric_finder.validation import (
    calculate_overlap,
    calculate_coverage,
    calculate_distance,
    validate_prediction
)


def test_validation():

    predicted_start = 3924600
    predicted_end = 3926500

    known_start = 3925597
    known_end = 3926500

    result = validate_prediction(
        predicted_start,
        predicted_end,
        known_start,
        known_end
    )

    print(result)

    assert result["overlap"] == 903
    assert result["distance"] == 0
