"""Tests for the bounded least-squares fitter (fit_least_squares).

This is the fast, bounded optimizer that makes large fitting sweeps feasible. The tests check that it
recovers identifiable parameters, reaches a fit comparable to the slower Nelder-Mead optimizer, and
honours its parameter bounds.
"""
import numpy as np

from mechanistic_inference import fit_least_squares, fit_mle, observed
from mechanistic_inference import models as M

T = np.linspace(0.0, 300.0, 60)
SIGMA = 0.03
TRUE3 = np.log([1.04, 0.35, 7.70, 0.02, 0.02])


def _setup():
    base = M.make_dimensional_simulator("M3")
    clean = base(TRUE3, T)
    scale = np.array([clean[:, 0].max(), clean[:, 1].max()])
    sim = lambda th, t: base(th, t) / scale  # noqa: E731
    return sim, clean, scale


def _rss(sim, theta, ch, data):
    return float(np.sum((observed(sim(theta, T), ch) - data) ** 2))


def test_least_squares_recovers_identifiable_params():
    """With both channels (identifiable), the fit recovers all parameters."""
    sim, clean, scale = _setup()
    rng = np.random.default_rng(0)
    data = (clean / scale)[:, [0, 1]].ravel() + rng.normal(scale=SIGMA, size=2 * T.size)
    mle = fit_least_squares(T, data, (0, 1), SIGMA, TRUE3 + 0.3, simulate=sim)
    assert np.max(np.abs(mle - TRUE3)) < 0.2


def test_least_squares_matches_nelder_mead_fit():
    """On the same data, the bounded least-squares fit reaches a residual sum of squares comparable to
    the Nelder-Mead optimizer (it is faster, not worse)."""
    sim, clean, scale = _setup()
    rng = np.random.default_rng(1)
    data = (clean / scale)[:, [1]].ravel() + rng.normal(scale=SIGMA, size=T.size)
    ls = fit_least_squares(T, data, (1,), SIGMA, TRUE3, simulate=sim)
    nm = fit_mle(T, data, (1,), SIGMA, TRUE3, simulate=sim, maxiter=400)
    rss_ls, rss_nm = _rss(sim, ls, (1,), data), _rss(sim, nm, (1,), data)
    assert rss_ls <= rss_nm + 0.05 * rss_nm + 1e-3


def test_multistart_m3_respects_nesting():
    """M1 is M3 with kappa_P = 0, so a multi-start M3 fit must not fit worse than M1 (residual sum of
    squares of M3 <= that of M1). A single-start M3 fit can violate this by getting stuck near kappa = 0;
    this guards against that bug returning in the discrimination experiments."""
    base3 = M.make_dimensional_simulator("M3")
    clean = base3(TRUE3, T)
    scale = np.array([clean[:, 0].max(), clean[:, 1].max()])
    sim3 = lambda th, t: base3(th, t) / scale  # noqa: E731
    base1 = M.make_dimensional_simulator("M1")
    sim1 = lambda th, t: base1(th, t) / scale  # noqa: E731
    rng = np.random.default_rng(0)
    data = (clean / scale)[:, [1]].ravel() + rng.normal(scale=SIGMA, size=T.size)

    mle1 = fit_least_squares(T, data, (1,), SIGMA, np.log([1.04, 0.35, 7.70, 0.02]), simulate=sim1)
    rss1 = _rss(sim1, mle1, (1,), data)
    best3 = np.inf
    for k in (1e-3, 5e-3, 2e-2, 8e-2, 3e-1):  # multi-start over kappa, including the M1-fit-based start
        m = fit_least_squares(T, data, (1,), SIGMA, np.append(mle1, np.log(k)), simulate=sim3)
        best3 = min(best3, _rss(sim3, m, (1,), data))
    assert best3 <= rss1 + 1e-3


def test_least_squares_respects_bounds():
    """The fit stays within theta0 +/- log_bound, which is what keeps it out of the slow stiff regions."""
    sim, clean, scale = _setup()
    rng = np.random.default_rng(2)
    data = (clean / scale)[:, [1]].ravel() + rng.normal(scale=SIGMA, size=T.size)
    mle = fit_least_squares(T, data, (1,), SIGMA, TRUE3, simulate=sim, log_bound=2.0)
    assert np.all(mle >= TRUE3 - 2.0 - 1e-6)
    assert np.all(mle <= TRUE3 + 2.0 + 1e-6)
