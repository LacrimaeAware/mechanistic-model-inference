"""Phase 3 extension: MCMC posterior for the regulated model M3 (no closed form, so the ODE is used).

The M1 posterior (recovery_and_posterior.py) showed the protein-only degeneracy as a ridge. This checks
the same in a regulated model and the kappa_P marginal (expected wide under protein-only, tight with
mRNA, matching the marginal practical identifiability of kappa).

Sampling note (learned the hard way): under protein-only the posterior is a long flat valley along the
non-identifiable direction, truncated by the prior. emcee's stretch move is affine-invariant, so
rotating coordinates does not help; the valley is non-elliptical, so it simply needs many steps to fill
(autocorrelation time ~200, so ~15000 steps for steps/tau > 50). The ODE makes each step expensive, so
the fix is to make each step cheap: a coarse time grid and a loose solver tolerance, which the posterior
shape is insensitive to. A built-in convergence check (steps / autocorrelation time) labels the result;
the metrics are only reported as final if it clears the threshold.

Run: python experiments/03_inference/mcmc_m3.py
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
from mechanistic_inference import neg_log_likelihood  # noqa: E402

RESULTS = os.path.join(os.path.dirname(__file__), "..", "..", "results")
os.makedirs(RESULTS, exist_ok=True)
T = np.linspace(0.0, 300.0, 20)  # coarse grid: makes each ODE solve cheap so many steps are affordable
SIGMA = 0.03
TRUE = np.log([1.04, 0.35, 7.70, 0.02, 0.02])  # [log k_m, d_m, k_p, d_p, kappa_P]
lines = []


def record(s):
    print(s); lines.append(s)


def main():
    try:
        import emcee
    except Exception as e:
        record(f"emcee unavailable ({e}); aborting.")
        return

    base = M.make_dimensional_simulator("M3", rtol=1e-3, atol=1e-6)  # loose tol; shape is insensitive
    clean = base(TRUE, T)
    scale = np.array([clean[:, 0].max(), clean[:, 1].max()])
    sim_norm = lambda th, t: base(th, t) / scale  # noqa: E731
    lo, hi = TRUE - 4.0, TRUE + 4.0

    def log_prob(theta, data, ch):
        if np.any(theta < lo) or np.any(theta > hi):
            return -np.inf
        val = -neg_log_likelihood(theta, T, data, ch, SIGMA, simulate=sim_norm)
        return val if np.isfinite(val) else -np.inf

    record("M3 MCMC posterior (emcee), protein-only vs mRNA+protein")
    record("=" * 64)
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2))
    # protein-only is the hard flat-ridge case and gets the longer chain
    settings = {"protein-only": (20000, 6000), "mRNA+protein": (8000, 2500)}
    for ax, (ch, lbl) in zip(axes, [((1,), "protein-only"), ((0, 1), "mRNA+protein")]):
        nsteps, burn = settings[lbl]
        rng = np.random.default_rng(0)
        data = (clean / scale)[:, list(ch)].ravel() + rng.normal(scale=SIGMA, size=len(ch) * T.size)
        nwalkers, ndim = 32, 5
        # independent tight ball; the many steps (not the init) fill the flat ridge. A collinear
        # ridge-seeded init makes the ensemble degenerate (large condition number), so do not do that.
        p0 = TRUE + 1e-2 * rng.normal(size=(nwalkers, ndim))
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args=(data, ch))
        sampler.run_mcmc(p0, nsteps, progress=False)
        s = sampler.get_chain(discard=burn, thin=10, flat=True)
        corr = np.corrcoef(s[:, 0], s[:, 2])[0, 1]
        ratio = np.std(s[:, 0] - s[:, 2]) / np.std(s[:, 0] + s[:, 2])
        acc = float(np.mean(sampler.acceptance_fraction))
        try:
            tau = float(np.mean(sampler.get_autocorr_time(tol=0)))
        except Exception:
            tau = float("nan")
        ess = (nsteps - burn) / tau if tau == tau else float("nan")
        conv = "converged" if ess > 50 else "UNDER-CONVERGED (preliminary)"
        record(f"  {lbl:12s}: corr(log k_m, log k_p)={corr:+.2f}  ridge ratio(diff/sum)={ratio:4.1f}  "
               f"kappa_P posterior sd={np.std(s[:, 4]):.3f}")
        record(f"               accept={acc:.2f}  autocorr~{tau:.0f}  steps/tau={ess:.0f}  -> {conv}")
        ax.scatter(s[:, 0] - TRUE[0], s[:, 2] - TRUE[2], s=4, alpha=0.2)
        ax.set_xlabel(r"$\log k_m$ offset"); ax.set_ylabel(r"$\log k_p$ offset")
        ax.set_title(f"M3 {lbl} (corr={corr:+.2f})")
        ax.set_xlim(-3, 3); ax.set_ylim(-3, 3)
    fig.suptitle("M3 posterior: protein-only ridge in (k_m, k_p), closed by observing mRNA")
    fig.tight_layout()
    f = os.path.join(RESULTS, "fig06_m3_posterior.png")
    fig.savefig(f, dpi=140)
    record(f"  wrote {os.path.relpath(f)}")
    out = os.path.join(RESULTS, "mcmc_m3_summary.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    record(f"  wrote {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
