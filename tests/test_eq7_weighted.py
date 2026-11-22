import numpy as np
import pandas as pd

from src.eq7_weighted import RecoveryPanel, fit_ridge, top_weights, training_contrasts


def counts(values):
    return pd.Series(values, index=pd.Index(
        [("C1", "V1"), ("C2", "V1"), ("C3", "V1"), ("C4", "V2")],
        tupleize_cols=False,
    )).loc[lambda s: s.gt(0)]


def test_recovery_uses_mouse_prevalence_and_full_library_denominator():
    train = RecoveryPanel({
        ("g1", "A"): counts([2, 4, 1, 0]),
        ("g1", "B"): counts([2, 0, 1, 0]),
        ("g5", "C"): counts([0, 0, 3, 4]),
    })
    test = RecoveryPanel({("g1", "D"): counts([2, 4, 1, 5])}, align_to=train)

    recovered = test.measure(train.reference(), [0])

    np.testing.assert_allclose(recovered[:, :, 0], [[200 / 3, 1e6 * (2 + 4 * 0.5) / 12]])


def test_outer_mouse_does_not_change_inner_response():
    inputs = {
        (group, str(i)): counts([i + 1, 3 - i, 0 if group == "g1" else 1, i + 1])
        for group in ("g1", "g5", "g6") for i in range(3)
    }
    first = RecoveryPanel(inputs)
    omitted = tuple(idx for idx, (group, mouse) in enumerate(first.mice) if mouse == "2")
    changed = {key: counts([0, 2, 20, 0]) if key[1] == "2" else row for key, row in inputs.items()}
    second = RecoveryPanel(changed)

    a = training_contrasts(first, omitted, {}).set_index("v_gene")
    b = training_contrasts(second, omitted, {}).set_index("v_gene")

    pd.testing.assert_frame_equal(a, b)


def test_signed_ridge_can_decrease_and_top_k_ties_are_fractional():
    frame = pd.DataFrame({"x": [1, 2, 3], "y": [3, 2, 1]})
    assert fit_ridge(frame, ("x",), "y", 1).coefficients["x"] < 0
    assert fit_ridge(frame, ("x",), "y", 1, nonnegative=True).coefficients["x"] == 0
    np.testing.assert_allclose(top_weights(np.array([4, 2, 2, 1]), 2), [1, 0.5, 0.5, 0])
