"""Experiment 06: discriminability map for transcriptional (M2) vs translational (M3) autoregulation.

Question (the one the source paper raises by comparing the two mechanisms, but does not answer): across
feedback strength and observation channel, how often can you correctly identify which mechanism generated
the data? Output is a selection rate per cell: the fraction of noise draws in which the true model wins
on AIC. ~50% means indistinguishable; ~100% means reliably distinguishable.

PROVISIONAL, first coarse pass, not independently re-verified. Lessons baked in: multi-start fitting
(single-start gets stuck near kappa = 0), live progress logging, rates over noise draws (never one fit).
Verification targets for this experiment are listed in docs/verification_debt.md.

Run: python -u experiments/06_discriminability_map/discriminability_map.py
"""
from __future__ import annotations

import os
import sys
import time
import warnings

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from mechanistic_inference import models as M  # noqa: E402
from mechanistic_inference import fit_least_squares, observed  # noqa: E402

RESULTS = os.path.join(os.path.dirname(__file__), "..", "..", "results")
os.makedirs(RESULTS, exist_ok=True)

T = np.linspace(0.0, 300.0, 60)
SIGMA = 0.03
N = 25                                  # noise draws per cell (a ~50% rate needs enough draws to pin)
MAX_NFEV = 200
RATES_BASE = [1.04, 0.35, 7.70, 0.02]   # k_m, d_m, k_p, d_p
KAPPAS = [0.005, 0.02, 0.05]            # 1/kappa = 200, 50, 20 molecules (weak, medium, strong feedback)
INV = {0.005: 200, 0.02: 50, 0.05: 20}
SCHEMES = [((1,), "protein-only"), ((0, 1), "mRNA+protein")]
KSTARTS = (1e-3, 5e-3, 2e-2, 8e-2, 3e-1)
lines = []


def log(s):
    print(s, flush=True); lines.append(s)


def truth_theta(model, kappa):
    return np.log(RATES_BASE + [kappa])  # M2 or M3 both take a single kappa as the 5th parameter


def fit_aic(candidate, data, ch, scale):
    """Multi-start (over kappa) AIC of the best fit for a regulated candidate model."""
    base = M.make_dimensional_simulator(candidate)
    sim = lambda th, t: base(th, t) / scale  # noqa: E731
    best = np.inf
    for k in KSTARTS:
        mle = fit_least_squares(T, data, ch, SIGMA, np.log(RATES_BASE + [k]), simulate=sim, max_nfev=MAX_NFEV)
        best = min(best, float(np.sum((observed(sim(mle, T), ch) - np.asarray(data)) ** 2)))
    return 2 * 5 + best / SIGMA ** 2


def select_rate(truth, kappa, ch, cell_id):
    other = "M3" if truth == "M2" else "M2"
    base = M.make_dimensional_simulator(truth)
    clean = base(truth_theta(truth, kappa), T)
    scale = np.array([clean[:, 0].max(), clean[:, 1].max()])
    wins = 0
    for d in range(N):
        rng = np.random.default_rng(cell_id * 1000 + d)   # reproducible integer seed (no hash())
        data = (clean / scale)[:, list(ch)].ravel() + rng.normal(scale=SIGMA, size=len(ch) * T.size)
        wins += int(fit_aic(truth, data, ch, scale) <= fit_aic(other, data, ch, scale))
    return wins / N


def main():
    log(f"Discriminability of transcriptional (M2) vs translational (M3) autoregulation.")
    log(f"Selection rate = fraction of {N} noise draws the true model wins on AIC (multi-start, "
        f"sigma={SIGMA}, grid={len(T)}). PROVISIONAL.")
    log("=" * 78)
    t0 = time.time()
    grid = {}  # (truth, kappa, scheme_label) -> rate
    cell_id = 0
    for truth in ("M2", "M3"):
        for kappa in KAPPAS:
            for ch, lbl in SCHEMES:
                cell_id += 1
                r = select_rate(truth, kappa, ch, cell_id)
                grid[(truth, kappa, lbl)] = r
                log(f"  truth {truth}, 1/kappa={INV[kappa]:>3} molec, {lbl:12s}: "
                    f"true model wins {r:4.0%}  [{time.time() - t0:.0f}s]")
        log("")
    log(f"done in {time.time() - t0:.0f}s")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, truth in zip(axes, ("M2", "M3")):
        Z = np.array([[grid[(truth, k, lbl)] for (_, lbl) in SCHEMES] for k in KAPPAS])
        im = ax.imshow(Z, vmin=0.5, vmax=1.0, cmap="viridis", aspect="auto")
        ax.set_xticks(range(len(SCHEMES))); ax.set_xticklabels([l for _, l in SCHEMES])
        ax.set_yticks(range(len(KAPPAS))); ax.set_yticklabels([f"1/k={INV[k]}" for k in KAPPAS])
        ax.set_title(f"truth = {truth} ({'transcriptional' if truth == 'M2' else 'translational'})")
        for i in range(len(KAPPAS)):
            for j in range(len(SCHEMES)):
                ax.text(j, i, f"{Z[i, j]:.0%}", ha="center", va="center",
                        color="w" if Z[i, j] < 0.8 else "k")
        fig.colorbar(im, ax=ax, label="true-model select rate")
    fig.suptitle("Discriminability of M2 vs M3: where, and from which channel (provisional)")
    fig.tight_layout()
    f = os.path.join(RESULTS, "fig10_discriminability_map.png")
    fig.savefig(f, dpi=140)
    log(f"wrote {os.path.relpath(f)}")
    out = os.path.join(RESULTS, "discriminability_map_summary.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    log(f"wrote {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
