"""Phase 2 methodology audit: check that the identifiability verdicts are real, not method artifacts.

Five checks, each targeting a way the structural/practical results could be fooling us:
  1. Spectral gap     - is the Fisher rank a clear call (big gap) or borderline?
  2. Step robustness  - is the near-zero eigenvalue real, or finite-difference noise? Does the null
                        direction match the analytic prediction [k_m +, k_p -]?
  3. Method agreement - does the profile-likelihood interval width agree with the Fisher standard
                        error for an identifiable parameter (two independent methods)?
  4. MLE recovery     - does fitting recover the truth where identifiable, and recover only the
                        k_m*k_p product (not the split) under protein-only?
  5. Noise sweep      - k_m should be non-identifiable at every noise level (structural); kappa should
                        cross from identifiable to not as noise grows (practical). Report the crossover.

Run: python experiments/02_identifiability/verify_methodology.py
"""
from __future__ import annotations

import os
import sys
import warnings

import numpy as np

# The optimizer probes extreme parameters at high noise, where the stiff solver issues convergence
# warnings; those regions are handled by the simulator's sentinel return, so silence the noise here.
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from mechanistic_inference import models as M  # noqa: E402
from mechanistic_inference import fisher_information, fit_mle, profile_likelihood, is_identifiable  # noqa: E402

RESULTS = os.path.join(os.path.dirname(__file__), "..", "..", "results")
os.makedirs(RESULTS, exist_ok=True)

TRUE = {
    "M1": np.log([1.04, 0.35, 7.70, 0.02]),
    "M3": np.log([1.04, 0.35, 7.70, 0.02, 0.02]),
}
NAMES3 = ["k_m", "d_m", "k_p", "d_p", "kappa_P"]
T = np.linspace(0.0, 300.0, 80)
lines = []


def record(s):
    print(s); lines.append(s)


def check_1_spectral_gap():
    record("CHECK 1: spectral gap (is the rank a clear call?)")
    for model, names in [("M1", ["k_m", "d_m", "k_p", "d_p"]), ("M3", NAMES3)]:
        sim = M.make_dimensional_simulator(model)
        for ch, lbl in [((1,), "protein-only"), ((0, 1), "mRNA+protein")]:
            ev = np.sort(np.clip(np.linalg.eigvalsh(fisher_information(TRUE[model], T, ch, 1.0, simulate=sim)), 0, None))
            rel = ev / ev.max()
            # gap = ratio between smallest "kept" and largest "dropped" eigenvalue
            kept = rel > 1e-9
            denom = rel[~kept].max() if (~kept).any() else 0.0
            gap = float("inf") if denom == 0.0 else rel[kept].min() / denom
            record(f"  {model} [{lbl:12s}] rel-eigs={np.array2string(rel, precision=1, formatter={'float': lambda x: f'{x:.1e}'})}")
            record(f"      smallest kept / largest dropped = {gap:.1e}  ({'clear gap' if gap > 1e3 else 'BORDERLINE'})")
    record("")


def check_2_step_robustness():
    record("CHECK 2: finite-difference step robustness and analytic null direction (M1 protein-only)")
    sim = M.make_dimensional_simulator("M1")
    analytic_null = np.array([1, 0, -1, 0]) / np.sqrt(2)  # log(k_m) - log(k_p) free, product fixed
    for eps in [1e-4, 1e-5, 1e-6]:
        F = fisher_information(TRUE["M1"], T, (1,), 1.0, simulate=sim, eps=eps)
        w, V = np.linalg.eigh(F)
        null = V[:, 0]
        align = abs(float(np.dot(null, analytic_null)))
        record(f"  eps={eps:.0e}: smallest eig (rel) = {w[0] / w[-1]:.1e}, |null . analytic| = {align:.4f}")
    record("  (alignment ~1.0 across step sizes means the null direction is the analytic k_m/k_p one,")
    record("   not numerical noise.)")
    record("")


def check_3_method_agreement():
    record("CHECK 3: profile interval width vs Fisher standard error (M3, mRNA+protein, identifiable case)")
    model = "M3"
    sim = M.make_dimensional_simulator(model)
    th = TRUE[model]
    clean = sim(th, T)
    scale = np.array([clean[:, 0].max(), clean[:, 1].max()])

    def sim_norm(theta, t):
        return M.make_dimensional_simulator(model)(theta, t) / scale

    sigma = 0.03
    ch = (0, 1)
    # Fisher standard errors (invert the full-rank FIM, in normalized units)
    F = fisher_information(th, T, ch, sigma, simulate=sim_norm)
    se = np.sqrt(np.diag(np.linalg.inv(F)))
    rng = np.random.default_rng(0)
    data = (clean / scale)[:, list(ch)].ravel()
    data = data + rng.normal(scale=sigma, size=data.size)
    mle = fit_mle(T, data, ch, sigma, th + 0.1, simulate=sim_norm, maxiter=800)
    for pidx in [0, 2, 4]:  # k_m, k_p, kappa_P
        grid, nll = profile_likelihood(pidx, mle, T, data, ch, sigma, simulate=sim_norm, span=2.0, n=21, maxiter=400)
        d = nll - nll.min()
        # half-width where the profile crosses 1.92 (one-parameter 95%), interpolated on the right side
        above = np.where((grid > mle[pidx]) & (d > 1.92))[0]
        half = (grid[above[0]] - mle[pidx]) if above.size else float("nan")
        record(f"  {NAMES3[pidx]:8s}: Fisher SE = {se[pidx]:.3f} (log units),  profile half-width(95%) ~ {half:.3f}")
    record("  (Fisher SE and 1.92-crossing half-width should be the same order; large disagreement would")
    record("   mean the quadratic Fisher approximation is misleading here.)")
    record("")


def check_4_mle_recovery():
    record("CHECK 4: MLE recovery (M3)")
    model = "M3"
    sim = M.make_dimensional_simulator(model)
    th = TRUE[model]
    clean = sim(th, T)
    scale = np.array([clean[:, 0].max(), clean[:, 1].max()])

    def sim_norm(theta, t):
        return M.make_dimensional_simulator(model)(theta, t) / scale

    sigma = 0.02
    rng = np.random.default_rng(1)
    # both channels: expect full recovery
    ch = (0, 1)
    data = (clean / scale)[:, list(ch)].ravel() + rng.normal(scale=sigma, size=T.size * 2)
    mle = fit_mle(T, data, ch, sigma, th + np.array([0.3, -0.2, 0.3, 0.1, -0.2]), simulate=sim_norm, maxiter=1500)
    record(f"  both channels: max |log param error| = {np.max(np.abs(mle - th)):.3f} (small = recovered)")
    record(f"    true   = {np.array2string(np.exp(th), precision=3)}")
    record(f"    fitted = {np.array2string(np.exp(mle), precision=3)}")
    # protein only: start offset ALONG the ridge (k_m up, k_p down); expect product preserved, split not pulled back
    ch = (1,)
    data = (clean / scale)[:, list(ch)].ravel() + rng.normal(scale=sigma, size=T.size)
    start = th + np.array([0.6, 0.0, -0.6, 0.0, 0.0])
    mle = fit_mle(T, data, ch, sigma, start, simulate=sim_norm, maxiter=1500)
    prod_err = abs((mle[0] + mle[2]) - (th[0] + th[2]))
    split_err = abs(mle[0] - th[0])
    record(f"  protein-only: |log(k_m*k_p) error| = {prod_err:.3f} (small = product recovered)")
    record(f"                |log k_m error|        = {split_err:.3f} (large = split NOT recovered -> ridge)")
    record(f"                kappa_P error          = {abs(mle[4] - th[4]):.3f} (small = kappa still recovered)")
    record("")


def check_5_noise_sweep():
    record("CHECK 5: noise sweep (k_m flat at all noise; kappa crosses identifiable -> not)")
    model = "M3"
    sim = M.make_dimensional_simulator(model)
    th = TRUE[model]
    clean = sim(th, T)
    scale = np.array([clean[:, 0].max(), clean[:, 1].max()])

    def sim_norm(theta, t):
        return M.make_dimensional_simulator(model)(theta, t) / scale

    record(f"  {'sigma':>7} | {'k_m (prot-only)':>16} | {'kappa_P (prot-only)':>20} | {'kappa_P (both)':>16}")
    for sigma in [0.005, 0.01, 0.03, 0.1, 0.3]:
        rng = np.random.default_rng(7)
        verdicts = {}
        for ch, key in [((1,), "km_p"), ((1,), "kap_p"), ((0, 1), "kap_b")]:
            data = (clean / scale)[:, list(ch)].ravel() + rng.normal(scale=sigma, size=len(ch) * T.size)
            mle = fit_mle(T, data, ch, sigma, th + 0.1, simulate=sim_norm, maxiter=500)
            pidx = 0 if key == "km_p" else 4
            _, nll = profile_likelihood(pidx, mle, T, data, ch, sigma, simulate=sim_norm, span=2.5, n=11, maxiter=300)
            verdicts[key] = "identifiable" if is_identifiable(nll) else "flat"
        record(f"  {sigma:>7.3f} | {verdicts['km_p']:>16} | {verdicts['kap_p']:>20} | {verdicts['kap_b']:>16}")
    record("  (k_m should read 'flat' at every sigma: non-identifiability is structural, not noise-driven.")
    record("   kappa_P should degrade to 'flat' as sigma grows: that is practical non-identifiability.)")
    record("")


def main():
    check_1_spectral_gap()
    check_2_step_robustness()
    check_3_method_agreement()
    check_4_mle_recovery()
    check_5_noise_sweep()
    out = os.path.join(RESULTS, "verify_methodology_summary.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    record(f"wrote {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
