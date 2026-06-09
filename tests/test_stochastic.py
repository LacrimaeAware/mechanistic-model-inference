"""Tests for the Gillespie SSA (src/mechanistic_inference/stochastic.py).

Seeds are fixed so the tests are deterministic. The CV separations between models are large relative
to ensemble noise, so the ordering assertions are robust; the exact CV values are not asserted (see
the propensity caveat in the module docstring).
"""
import numpy as np

from mechanistic_inference import stochastic as S

BASE = dict(alpha_M=1.04, beta_M=0.35, beta_P=0.02, n_M=2.0, n_P=2.0, S_a=0.0, M0=3, P0=400)
SAMPLE_TIMES = [250, 300, 350, 400, 450, 500, 550, 600]
KAP = 1.0 / 300.0


def test_ssa_mean_matches_ode_steady_state():
    """The SSA mean protein for the simplistic model sits near its deterministic steady state."""
    r = S.ensemble_cv(n_runs=120, sample_times=SAMPLE_TIMES, seed=0, alpha_P=2.94, **BASE)
    p_ode = 2.94 * (1.04 / 0.35) / 0.02  # alpha_P M* / beta_P, with M* = alpha_M / beta_M
    assert abs(r["mean_P"] - p_ode) / p_ode < 0.05


def test_protein_noise_ordering():
    """Protein CV orders as transcriptional > simplistic > translational (the paper's claim)."""
    m1 = S.ensemble_cv(n_runs=120, sample_times=SAMPLE_TIMES, seed=1, alpha_P=2.94, **BASE)
    m2 = S.ensemble_cv(n_runs=120, sample_times=SAMPLE_TIMES, seed=2, alpha_P=7.70, kappa_M=KAP, **BASE)
    m3 = S.ensemble_cv(n_runs=120, sample_times=SAMPLE_TIMES, seed=3, alpha_P=7.70, kappa_P=KAP, **BASE)
    assert m2["cv_P"] > m1["cv_P"] > m3["cv_P"]


def test_translational_does_not_change_mrna_noise():
    """M3 leaves mRNA unregulated, so its mRNA CV is close to the simplistic model's, while M2
    (transcriptional) raises mRNA noise."""
    m1 = S.ensemble_cv(n_runs=120, sample_times=SAMPLE_TIMES, seed=1, alpha_P=2.94, **BASE)
    m2 = S.ensemble_cv(n_runs=120, sample_times=SAMPLE_TIMES, seed=2, alpha_P=7.70, kappa_M=KAP, **BASE)
    m3 = S.ensemble_cv(n_runs=120, sample_times=SAMPLE_TIMES, seed=3, alpha_P=7.70, kappa_P=KAP, **BASE)
    assert m2["cv_M"] > m1["cv_M"]
    assert abs(m3["cv_M"] - m1["cv_M"]) < 0.1
