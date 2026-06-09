"""Correctness tests for the three negative-autoregulation models (src/mechanistic_inference/models.py).

Each test encodes a fact stated in Ryzowicz and Yildirim (2023) or a structural property that the
later identifiability work depends on. Parameter values used for the oscillation tests were checked
against the model's own discriminant (Eq. 9); see experiments/01.
"""
import numpy as np

from mechanistic_inference import models as M


def test_simplistic_matches_analytic():
    """M1 numerical integration matches the closed form (Eq. 8) across signal and parameter regimes."""
    t = np.linspace(0, 50, 800)
    for theta, S_a, m0, p0 in [(0.0571, 0.0, 0.0, 0.0), (0.6, 3.0, 0.2, 0.5), (0.2, 1.0, 1.0, 2.0)]:
        num = M.simulate_simplistic(t, theta=theta, S_a=S_a, m0=m0, p0=p0)
        ana = M.simplistic_analytic(t, theta=theta, S_a=S_a, m0=m0, p0=p0)
        assert np.max(np.abs(num - ana)) < 1e-5


def test_simplistic_steady_state():
    """M1 settles at (1 + S_a, (1 + S_a)/theta)."""
    theta, S_a = 0.0571, 3.0
    m_star, p_star = M.steady_state(theta=theta, S_a=S_a)
    assert np.isclose(m_star, 1.0 + S_a, rtol=1e-6)
    assert np.isclose(p_star, (1.0 + S_a) / theta, rtol=1e-6)


def test_translational_mrna_equals_simplistic_mrna():
    """M3 leaves mRNA unregulated, so its mRNA trajectory is identical to M1's.

    Consequence for identifiability: the mRNA channel cannot distinguish M1 from M3; only protein can.
    """
    t = np.linspace(0, 40, 400)
    m1 = M.simulate_simplistic(t, theta=0.0571, S_a=2.0)
    m3 = M.simulate_translational_nar(t, theta=0.0571, kappa_tilde_P=2.0, n_P=2, S_a=2.0)
    assert np.allclose(m3[:, 0], m1[:, 0], atol=1e-7)


def test_regulation_lowers_protein_steady_state():
    """Both autoregulation models sit at a lower protein steady state than the simplistic model."""
    theta, S_a = 0.0571, 2.0
    p1 = M.steady_state(theta=theta, S_a=S_a)[1]
    p2 = M.steady_state(theta=theta, kappa_tilde_M=2.0, n_M=2, S_a=S_a)[1]
    p3 = M.steady_state(theta=theta, kappa_tilde_P=2.0, n_P=2, S_a=S_a)[1]
    assert p2 < p1 and p3 < p1


def test_transcriptional_and_translational_coincide_when_symmetric():
    """M2 and M3 share the same protein steady state when kappa_tilde_M = kappa_tilde_P and n_M = n_P.

    This is the model-discrimination degeneracy: at symmetric feedback the two mechanisms are
    indistinguishable in the protein steady-state observable despite different dynamics.
    """
    p2 = M.steady_state(theta=0.3, kappa_tilde_M=1.5, n_M=2, S_a=1.0)[1]
    p3 = M.steady_state(theta=0.3, kappa_tilde_P=1.5, n_P=2, S_a=1.0)[1]
    assert np.isclose(p2, p3, rtol=1e-8)


def test_only_transcriptional_oscillates():
    """At theta=0.6 with strong inhibition (1/kappa_tilde_M = 0.5), M2 has complex eigenvalues
    (Delta < 0) and its protein overshoots; M1 and M3 never overshoot."""
    theta = 0.6
    assert M.m2_oscillation_discriminant(theta=theta, kappa_tilde_M=2.0, n_M=8) < 0
    t = np.linspace(0, 60, 3000)
    p2 = M.simulate_transcriptional_nar(t, theta=theta, kappa_tilde_M=2.0, n_M=8, S_a=0.0)[:, 1]
    p1 = M.simulate_simplistic(t, theta=theta, S_a=0.0)[:, 1]
    p3 = M.simulate_translational_nar(t, theta=theta, kappa_tilde_P=2.0, n_P=8, S_a=0.0)[:, 1]
    assert p2.max() - p2[-1] > 1e-2          # M2 overshoots (damped oscillation)
    assert p1.max() - p1[-1] < 1e-6          # M1 monotone
    assert p3.max() - p3[-1] < 1e-6          # M3 monotone


def test_weak_inhibition_is_monotone():
    """At theta=0.6 with weak inhibition (1/kappa_tilde_M = 6), M2 has real eigenvalues (Delta >= 0)
    and does not overshoot."""
    theta = 0.6
    assert M.m2_oscillation_discriminant(theta=theta, kappa_tilde_M=1 / 6, n_M=8) > 0
    t = np.linspace(0, 60, 3000)
    p2 = M.simulate_transcriptional_nar(t, theta=theta, kappa_tilde_M=1 / 6, n_M=8, S_a=0.0)[:, 1]
    assert p2.max() - p2[-1] < 1e-6


def test_autoregulation_responds_faster_than_simplistic():
    """Both autoregulation models reach half their step response sooner than the simplistic model."""
    theta, S_a = 0.0571, 3.0
    kt = (1.04 * 7.70 / 0.35 ** 2) * (1.0 / 50.0)  # 1/kappa = 50 molecules
    t = np.linspace(0, 400, 4000)

    def half_time(p):
        target = p[0] + 0.5 * (p[-1] - p[0])
        idx = np.where((p - target) * np.sign(p[-1] - p[0]) >= 0)[0]
        return t[idx[0]]

    s1 = M.steady_state(theta=theta, S_a=0.0)
    s2 = M.steady_state(theta=theta, kappa_tilde_M=kt, n_M=2, S_a=0.0)
    s3 = M.steady_state(theta=theta, kappa_tilde_P=kt, n_P=2, S_a=0.0)
    th1 = half_time(M.simulate_simplistic(t, theta=theta, S_a=S_a, m0=s1[0], p0=s1[1])[:, 1])
    th2 = half_time(M.simulate_transcriptional_nar(t, theta=theta, kappa_tilde_M=kt, n_M=2, S_a=S_a, m0=s2[0], p0=s2[1])[:, 1])
    th3 = half_time(M.simulate_translational_nar(t, theta=theta, kappa_tilde_P=kt, n_P=2, S_a=S_a, m0=s3[0], p0=s3[1])[:, 1])
    assert th2 < th1 and th3 < th1
