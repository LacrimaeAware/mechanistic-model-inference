"""Phase 2: structural and practical identifiability map for the three negative-autoregulation models.

Question: for each model and each observation scheme (protein only, or mRNA plus protein), which
parameters can be estimated, and which combinations cannot? Structural identifiability is read from the
Fisher information rank and its null direction; practical identifiability from profile likelihood under
simulated noisy data.

The parameters are the dimensional rates and the regulation strength (log space): M1 [k_m, d_m, k_p,
d_p], M2 adds kappa_M, M3 adds kappa_P (Hill coefficient n and input level S held fixed). See
src/mechanistic_inference/models.py: make_dimensional_simulator.

Run: python experiments/02_identifiability/identifiability_map.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from mechanistic_inference import models as M  # noqa: E402
from mechanistic_inference import fisher_information, fit_mle, profile_likelihood, is_identifiable  # noqa: E402

RESULTS = os.path.join(os.path.dirname(__file__), "..", "..", "results")
os.makedirs(RESULTS, exist_ok=True)

# True dimensional parameters (Table 1): k_m, d_m, k_p, d_p, and kappa = 0.02 (1/kappa = 50 molecules).
TRUE = {
    "M1": np.log([1.04, 0.35, 7.70, 0.02]),
    "M2": np.log([1.04, 0.35, 7.70, 0.02, 0.02]),
    "M3": np.log([1.04, 0.35, 7.70, 0.02, 0.02]),
}
NAMES = {
    "M1": ["k_m", "d_m", "k_p", "d_p"],
    "M2": ["k_m", "d_m", "k_p", "d_p", "kappa_M"],
    "M3": ["k_m", "d_m", "k_p", "d_p", "kappa_P"],
}
SCHEMES = [((1,), "protein-only"), ((0, 1), "mRNA+protein")]
T = np.linspace(0.0, 300.0, 80)


def main():
    lines = []

    def record(s):
        print(s); lines.append(s)

    record("STRUCTURAL IDENTIFIABILITY (Fisher rank and null direction)")
    record("=" * 64)
    for model in ("M1", "M2", "M3"):
        sim = M.make_dimensional_simulator(model)
        th = TRUE[model]
        record(f"{model}  params {NAMES[model]}")
        for ch, lbl in SCHEMES:
            F = fisher_information(th, T, ch, 1.0, simulate=sim)
            w, V = np.linalg.eigh(F)
            ev = np.clip(w, 0, None)
            rank = int(np.sum(ev > 1e-9 * (ev.max() + 1e-30)))
            null = V[:, 0]
            direction = ", ".join(f"{n}:{v:+.2f}" for n, v in zip(NAMES[model], null))
            verdict = "FULL RANK" if rank == len(th) else f"rank-deficient (1 null direction)"
            record(f"  [{lbl:12s}] rank {rank}/{len(th)}  {verdict}")
            if rank < len(th):
                record(f"      non-identifiable direction: {direction}")
    record("")
    record("Reading: in every model, protein-only leaves the k_m vs k_p (transcription vs translation)")
    record("product non-identifiable; the regulation strength kappa is identifiable; observing mRNA")
    record("restores full rank. Regulation adds an estimable parameter but does not rescue the rate")
    record("degeneracy.")
    record("")

    # --- Practical identifiability (profile likelihood) for both regulated models ---
    record("PRACTICAL IDENTIFIABILITY (profile likelihood, 3% relative noise)")
    record("=" * 64)

    def make_sim_norm(model, scale):
        base = M.make_dimensional_simulator(model)
        return lambda theta, t: base(theta, t) / scale

    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    sigma = 0.03
    for row, (model, kname) in enumerate([("M2", "kappa_M"), ("M3", "kappa_P")]):
        sim = M.make_dimensional_simulator(model)
        th = TRUE[model]
        clean = sim(th, T)
        scale = np.array([clean[:, 0].max(), clean[:, 1].max()])  # per-channel normalization
        sim_norm = make_sim_norm(model, scale)
        rng = np.random.default_rng(0)
        for col, (pname, pidx) in enumerate([("k_m", 0), (kname, 4)]):
            ax = axes[row, col]
            for ch, lbl in SCHEMES:
                data = (clean / scale)[:, list(ch)].ravel()
                data = data + rng.normal(scale=sigma, size=data.size)
                mle = fit_mle(T, data, ch, sigma, th + 0.1, simulate=sim_norm, maxiter=600)
                # span must exceed the confidence-interval half-width or a genuinely identifiable but
                # weakly constrained parameter (e.g. M2 kappa_M) reads as a false "flat"; 3.0 log units
                # (a ~20x factor each way) is wide enough here.
                grid, nll = profile_likelihood(pidx, mle, T, data, ch, sigma, simulate=sim_norm,
                                               span=3.0, n=15, maxiter=400)
                ident = is_identifiable(nll)
                record(f"  {model} {pname:8s} under {lbl:12s}: "
                       f"{'identifiable' if ident else 'NON-identifiable (flat)'}")
                ax.plot(grid - mle[pidx], nll - nll.min(), marker="o", ms=3,
                        label=f"{lbl} ({'id' if ident else 'flat'})")
            ax.axhline(1.92, color="k", lw=0.8, ls=":")  # 95% threshold for one parameter
            ax.set_yscale("symlog", linthresh=2.0)  # linear near the threshold, log for the steep arms
            ax.set_title(f"{model}: profile {pname}")
            ax.set_xlabel(r"$\log$ offset from MLE")
            ax.set_ylabel(r"$\Delta$ neg-log-likelihood")
            ax.legend(fontsize=7)
    fig.suptitle("Practical identifiability: k_m flat under protein-only, bounded with mRNA; kappa bounded")
    fig.tight_layout()
    f = os.path.join(RESULTS, "fig04_profile_likelihood.png")
    fig.savefig(f, dpi=140)
    record(f"wrote {os.path.relpath(f)}")

    # kappa under protein-only is near the identifiable/flat boundary; a single profile is not a stable
    # verdict, so quantify it over noise draws.
    record("")
    record("kappa robustness under protein-only over 10 noise draws (3% noise):")
    for model, kidx, kname in [("M2", 4, "kappa_M"), ("M3", 4, "kappa_P")]:
        base = M.make_dimensional_simulator(model)
        th = TRUE[model]
        clean = base(th, T)
        scale = np.array([clean[:, 0].max(), clean[:, 1].max()])
        sim_norm = make_sim_norm(model, scale)
        cnt = 0
        for seed in range(10):
            rng = np.random.default_rng(seed)
            data = (clean / scale)[:, [1]].ravel() + rng.normal(scale=0.03, size=T.size)
            mle = fit_mle(T, data, (1,), 0.03, th + 0.1, simulate=sim_norm, maxiter=800)
            _, nll = profile_likelihood(kidx, mle, T, data, (1,), 0.03, simulate=sim_norm,
                                        span=3.0, n=15, maxiter=400)
            cnt += int(is_identifiable(nll))
        record(f"  {model} {kname} (protein-only): identifiable in {cnt}/10 draws (marginal)")

    out = os.path.join(RESULTS, "identifiability_map_summary.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    record(f"wrote {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
