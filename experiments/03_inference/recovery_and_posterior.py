"""Phase 3: inference under noise. Two parts.

A. Recovery vs noise (M3): repeated maximum-likelihood fits to simulated noisy data, reporting how well
   each parameter is recovered as a function of noise, under protein-only vs mRNA+protein observation.
   The structural prediction is that protein-only recovers the k_m*k_p product but not the split.

B. Posterior geometry (M1, MCMC): sample the Bayesian posterior with emcee and plot k_m vs k_p. The
   protein-only non-identifiability should appear as a ridge (a 1-D valley along log k_m + log k_p =
   const); observing mRNA should close it into a bounded blob. M1 uses its closed-form solution so the
   MCMC is fast and exact.

Run: python experiments/03_inference/recovery_and_posterior.py
"""
from __future__ import annotations

import os
import sys
import warnings

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from mechanistic_inference import models as M  # noqa: E402
from mechanistic_inference import fit_mle, neg_log_likelihood, observed  # noqa: E402

RESULTS = os.path.join(os.path.dirname(__file__), "..", "..", "results")
os.makedirs(RESULTS, exist_ok=True)
T = np.linspace(0.0, 300.0, 80)
lines = []


def record(s):
    print(s); lines.append(s)


def cascade_analytic(theta, t):
    """Closed-form M1 cascade (dimensional), theta = log[k_m, d_m, k_p, d_p], zero initial conditions."""
    k_m, d_m, k_p, d_p = np.exp(np.asarray(theta, dtype=float))
    t = np.asarray(t, dtype=float)
    m = (k_m / d_m) * (1.0 - np.exp(-d_m * t))
    amp = k_p * k_m / d_m
    term1 = (1.0 - np.exp(-d_p * t)) / d_p
    diff = d_p - d_m
    term2 = t * np.exp(-d_m * t) if abs(diff) < 1e-8 else (np.exp(-d_m * t) - np.exp(-d_p * t)) / diff
    p = amp * (term1 - term2)
    return np.column_stack([m, p])


def part_a_recovery():
    record("PART A: recovery vs noise (M3, repeated fits)")
    record("=" * 64)
    model = "M3"
    true = np.log([1.04, 0.35, 7.70, 0.02, 0.02])
    names = ["k_m", "d_m", "k_p", "d_p", "kappa_P"]
    base = M.make_dimensional_simulator(model)
    clean = base(true, T)
    scale = np.array([clean[:, 0].max(), clean[:, 1].max()])
    sim_norm = lambda th, t: base(th, t) / scale  # noqa: E731

    for ch, lbl in [((1,), "protein-only"), ((0, 1), "mRNA+protein")]:
        record(f"  scheme: {lbl}")
        for sigma in [0.01, 0.05]:
            errs, prod_errs, split_errs = [], [], []
            for r in range(12):
                rng = np.random.default_rng(100 + r)
                data = (clean / scale)[:, list(ch)].ravel() + rng.normal(scale=sigma, size=len(ch) * T.size)
                start = true + rng.normal(scale=0.2, size=true.size)
                mle = fit_mle(T, data, ch, sigma, start, simulate=sim_norm, maxiter=1200)
                errs.append(np.abs(mle - true))
                prod_errs.append(abs((mle[0] + mle[2]) - (true[0] + true[2])))
                split_errs.append(abs(mle[0] - true[0]))
            med = np.median(np.array(errs), axis=0)
            record(f"    sigma={sigma:.2f}: median |log error| per param = "
                   + ", ".join(f"{n}:{e:.3f}" for n, e in zip(names, med)))
            record(f"             k_m*k_p product err (median) = {np.median(prod_errs):.3f}; "
                   f"k_m split err (median) = {np.median(split_errs):.3f}")
    record("")


def part_b_posterior():
    record("PART B: posterior geometry (M1, MCMC)")
    record("=" * 64)
    try:
        import emcee
    except Exception as e:
        record(f"  emcee unavailable ({e}); skipping MCMC.")
        return

    true = np.log([1.04, 0.35, 7.70, 0.02])
    scale = np.array([1.04 / 0.35, 7.70 * 1.04 / (0.35 * 0.02)])  # steady-state [m, p]
    sim_norm = lambda th, t: cascade_analytic(th, t) / scale  # noqa: E731
    sigma = 0.03
    lo, hi = true - 4.0, true + 4.0

    def log_prob(theta, data, ch):
        if np.any(theta < lo) or np.any(theta > hi):
            return -np.inf
        return -neg_log_likelihood(theta, T, data, ch, sigma, simulate=sim_norm)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2))
    for ax, (ch, lbl) in zip(axes, [((1,), "protein-only"), ((0, 1), "mRNA+protein")]):
        rng = np.random.default_rng(0)
        clean = sim_norm(true, T)
        data = clean[:, list(ch)].ravel() + rng.normal(scale=sigma, size=len(ch) * T.size)
        nwalkers, ndim, nsteps, burn = 32, 4, 16000, 5000
        p0 = true + 1e-2 * rng.normal(size=(nwalkers, ndim))
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args=(data, ch))
        sampler.run_mcmc(p0, nsteps, progress=False)
        s = sampler.get_chain(discard=burn, thin=10, flat=True)
        corr = np.corrcoef(s[:, 0], s[:, 2])[0, 1]
        spread_sum = np.std(s[:, 0] + s[:, 2])   # along the product direction
        spread_diff = np.std(s[:, 0] - s[:, 2])  # along the ridge direction
        try:
            tau = float(np.mean(sampler.get_autocorr_time(tol=0)))
        except Exception:
            tau = float("nan")
        ess = (nsteps - burn) / tau if tau == tau else float("nan")
        conv = "converged" if ess > 50 else "under-converged"
        record(f"  {lbl:12s}: corr(log k_m, log k_p) = {corr:+.2f}; "
               f"std(sum)={spread_sum:.3f}, std(diff)={spread_diff:.3f}, "
               f"ratio diff/sum={spread_diff / spread_sum:.1f}; steps/tau={ess:.0f} ({conv})")
        ax.scatter(s[:, 0] - true[0], s[:, 2] - true[2], s=4, alpha=0.2)
        ax.set_xlabel(r"$\log k_m$ offset"); ax.set_ylabel(r"$\log k_p$ offset")
        ax.set_title(f"{lbl}  (corr={corr:+.2f})")
        ax.set_xlim(-3, 3); ax.set_ylim(-3, 3)
    fig.suptitle("M1 posterior: protein-only is a ridge along log k_m + log k_p = const; mRNA closes it")
    fig.tight_layout()
    f = os.path.join(RESULTS, "fig05_posterior_ridge.png")
    fig.savefig(f, dpi=140)
    record(f"  wrote {os.path.relpath(f)}")
    record("  (A long ridge -> ratio diff/sum >> 1 and corr ~ -1 means the product is constrained but")
    record("   the split is not. A round blob -> ratio ~ 1 and corr ~ 0 means both are constrained.)")
    record("")


def main():
    part_a_recovery()
    part_b_posterior()
    out = os.path.join(RESULTS, "inference_summary.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    record(f"wrote {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
