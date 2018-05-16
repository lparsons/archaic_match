from archaic_match.__main__ import calculate_thresholds


def test_calculate_thresholds_0():
    (lower, upper) = calculate_thresholds(0, 10)
    assert lower == 10
    assert upper == 10


def test_calculate_thresholds_1():
    (lower, upper) = calculate_thresholds(1, 10)
    assert lower == 9
    assert upper == 11


def test_calculate_thresholds_0_1():
    (lower, upper) = calculate_thresholds(0.1, 9)
    assert lower == 8.1
    assert upper == 9.9
