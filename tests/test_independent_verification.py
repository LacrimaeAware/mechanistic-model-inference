"""Independent re-derivations of the project's load-bearing findings, by a different route.

These are the "triple-check" tests: each confirms a headline claim with a method independent of the code
that produced the original result, so the conclusions do not rest on a single implementation.

- test_symbolic_protein_only_degeneracy: a symbolic (sympy) proof, using nothing from this project's
  tooling, that protein-only data sees k_m and k_p only through their product, and that mRNA breaks it.
- test_scaling_invariance_regulated_models: a high-precision numeric confirmation that the same
  degeneracy holds for the regulated (nonlinear) models M2 and M3.
- test_fisher_null_matches_proof: the numerical Fisher null direction equals the analytic one.
- test_distinguishability_multistart_is_global: the mechanism-mimic residual does not improve with more
  optimizer starts, so the distinguishability result is the global optimum, not an under-fit.
"""
import numpy as np

from mechanistic_inference import fisher_information, fit_least_squares, observed
from mechanistic_inference import models as M

T = np.linspace(0.0, 300.0, 200)
RATES = [1.04, 0.35, 7.70, 0.02]


def test_symbolic_protein_only_degeneracy():
    """Algebraic proof (independent of this project's numerics): in the mRNA->protein cascade the protein
    sees k_m and k_p only through their product, and observing mRNA carries information about k_m but not
    k_p, so it breaks the degeneracy."""
    import sympy as sp

    t, k_m, d_m, k_p, d_p = sp.symbols("t k_m d_m k_p d_p", positive=True)
    Mf, Pf = sp.Function("M"), sp.Function("P")
    m_sol = sp.dsolve(sp.Eq(Mf(t).diff(t), k_m - d_m * Mf(t)), Mf(t), ics={Mf(0): 0}).rhs
    p_sol = sp.dsolve(sp.Eq(Pf(t).diff(t), k_p * m_sol - d_p * Pf(t)), Pf(t), ics={Pf(0): 0}).rhs

    # protein sensitivity to log k_m equals that to log k_p => exactly one non-identifiable direction
    assert sp.simplify(k_m * sp.diff(p_sol, k_m) - k_p * sp.diff(p_sol, k_p)) == 0
    # and that sensitivity equals P itself => P is proportional to k_m (hence to the product k_m*k_p)
    assert sp.simplify(k_m * sp.diff(p_sol, k_m) - p_sol) == 0
    # mRNA carries no information about k_p, but does about k_m => observing mRNA breaks the tie
    assert sp.simplify(k_p * sp.diff(m_sol, k_p)) == 0
    assert sp.simplify(k_m * sp.diff(m_sol, k_m)) != 0


def test_scaling_invariance_regulated_models():
    """The protein trajectory of each regulated model is invariant under k_m -> c*k_m, k_p -> k_p/c
    (high-precision), confirming the product degeneracy holds for the nonlinear M2 and M3, not only M1."""
    c = 2.3
    for model in ("M2", "M3"):
        base = M.make_dimensional_simulator(model, rtol=1e-10, atol=1e-12)
        th0 = np.log(RATES + [0.02])
        th1 = th0.copy()
        th1[0] += np.log(c)
        th1[2] -= np.log(c)
        p0, p1 = base(th0, T)[:, 1], base(th1, T)[:, 1]
        assert np.max(np.abs(p0 - p1)) / p0.max() < 1e-6


def test_fisher_null_matches_proof():
    """The numerical Fisher null direction (protein-only) equals the analytic k_m-vs-k_p direction."""
    sim = M.make_dimensional_simulator("M1")
    F = fisher_information(np.log(RATES), T, (1,), 1.0, simulate=sim)
    _, V = np.linalg.eigh(F)
    assert abs(float(V[:, 0] @ (np.array([1, 0, -1, 0]) / np.sqrt(2)))) > 0.999


def test_distinguishability_multistart_is_global():
    """The best M2-mimics-M3 protein residual does not improve when the multi-start is doubled, so the
    distinguishability number is the global optimum (not an under-fit), and it is tiny (<1% of peak)."""
    c3 = M.make_dimensional_simulator("M3")(np.log(RATES + [0.02]), T)
    scale = np.array([c3[:, 0].max(), c3[:, 1].max()])
    data = observed(c3 / scale, (1,))
    base = M.make_dimensional_simulator("M2")
    sim = lambda th, t: base(th, t) / scale  # noqa: E731

    def best(ks):
        out = np.inf
        for k in ks:
            mle = fit_least_squares(T, data, (1,), 1.0, np.log(RATES + [k]), simulate=sim, max_nfev=600)
            out = min(out, float(np.sqrt(np.mean((observed(sim(mle, T), (1,)) - data) ** 2))))
        return out

    r5 = best((1e-3, 5e-3, 2e-2, 8e-2, 3e-1))
    r11 = best((5e-4, 1e-3, 3e-3, 5e-3, 1e-2, 2e-2, 4e-2, 8e-2, 1.5e-1, 3e-1, 5e-1))
    assert abs(r5 - r11) < 1e-4        # more starts does not help => global optimum reached
    assert r5 < 0.01                   # < 1% of peak => indistinguishable from protein
