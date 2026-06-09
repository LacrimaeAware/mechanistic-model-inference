"""Phase 2 structural-identifiability tests for the three negative-autoregulation models.

These check the Fisher-rank verdicts of the identifiability map and the analytic reason behind the
non-identifiable direction: the protein trajectory is invariant under the scaling k_m -> k_m c,
k_p -> k_p / c, so protein-only observation cannot separate transcription rate from translation rate,
in the regulated models as well as the simplistic one.
"""
import numpy as np

from mechanistic_inference import identifiability_report
from mechanistic_inference import models as M

TRUE = {
    "M1": np.log([1.04, 0.35, 7.70, 0.02]),
    "M2": np.log([1.04, 0.35, 7.70, 0.02, 0.02]),
    "M3": np.log([1.04, 0.35, 7.70, 0.02, 0.02]),
}
T = np.linspace(0.0, 300.0, 60)


def _rank(model, channels):
    sim = M.make_dimensional_simulator(model)
    rep = identifiability_report(TRUE[model], T, channels=channels, sigma=1.0, simulate=sim)
    return rep["rank"], rep["n_params"]


def test_m1_protein_only_rank_deficient():
    rank, n = _rank("M1", (1,))
    assert (rank, n) == (3, 4)  # k_m vs k_p product non-identifiable


def test_m1_both_channels_full_rank():
    assert _rank("M1", (0, 1)) == (4, 4)


def test_m2_protein_only_one_null_direction():
    rank, n = _rank("M2", (1,))
    assert (rank, n) == (4, 5)  # kappa_M identifiable, but k_m vs k_p still degenerate


def test_m2_both_channels_full_rank():
    assert _rank("M2", (0, 1)) == (5, 5)


def test_m3_protein_only_one_null_direction():
    rank, n = _rank("M3", (1,))
    assert (rank, n) == (4, 5)  # kappa_P identifiable, but k_m vs k_p still degenerate


def test_m3_both_channels_full_rank():
    assert _rank("M3", (0, 1)) == (5, 5)


def test_rate_product_degeneracy_leaves_protein_invariant():
    """Scaling k_m by c and k_p by 1/c leaves the protein trajectory identical while changing mRNA.

    This is the analytic reason protein-only observation cannot identify k_m and k_p separately, and it
    holds for both regulated models (the regulation strength is untouched by the scaling).
    """
    t = np.linspace(0.0, 300.0, 120)
    c = 1.7
    for model, kappa in [("M2", 0.02), ("M3", 0.02)]:
        sim = M.make_dimensional_simulator(model)
        base = np.log([1.04, 0.35, 7.70, 0.02, kappa])
        scaled = base.copy()
        scaled[0] += np.log(c)   # k_m * c
        scaled[2] -= np.log(c)   # k_p / c
        a = sim(base, t)
        b = sim(scaled, t)
        assert np.allclose(a[:, 1], b[:, 1], rtol=1e-4, atol=1e-4)      # protein invariant
        assert not np.allclose(a[:, 0], b[:, 0], rtol=1e-2)            # mRNA changes (scaled by c)
