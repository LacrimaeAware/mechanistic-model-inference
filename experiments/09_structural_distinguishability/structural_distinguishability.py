"""Experiment 09: structural distinguishability of transcriptional (M2) vs translational (M3) autoregulation.

The discrimination question done the RIGHT way: deterministically, with no noise draws, so the answer is
stable (the lesson from the failed AIC-rate experiments). For each feedback strength and observation
channel, fit each model to the OTHER model's noise-free output (best-effort mimicry, multi-start) and
measure the leftover residual as a percent of the signal. That residual is the structural distance
between the mechanisms: if it is below the measurement noise, the mechanisms are indistinguishable; if
above, they are distinguishable. This is the open question the source paper never addressed.

Run: python experiments/09_structural_distinguishability/structural_distinguishability.py
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
from mechanistic_inference import fit_least_squares, observed  # noqa: E402

RESULTS = os.path.join(os.path.dirname(__file__), "..", "..", "results")
os.makedirs(RESULTS, exist_ok=True)
T = np.linspace(0.0, 300.0, 120)
RATES = [1.04, 0.35, 7.70, 0.02]                    # k_m, d_m, k_p, d_p
KAPPAS = [0.005, 0.02, 0.05, 0.1, 0.2]              # 1/kappa = 200, 50, 20, 10, 5 molecules
INV = {0.005: 200, 0.02: 50, 0.05: 20, 0.1: 10, 0.2: 5}
KSTARTS = (1e-3, 5e-3, 2e-2, 8e-2, 3e-1)
lines = []


def log(s):
    print(s); lines.append(s)


def mimic_residual(fit_model, target, channels):
    """Best-effort (multi-start, noise-free) fit of fit_model to target's channels. Returns RMS residual
    as a percent of the target's peak (per-channel normalized so each channel contributes on its scale)."""
    scale = np.array([target[:, 0].max(), target[:, 1].max()])
    data = observed(target / scale, channels)
    base = M.make_dimensional_simulator(fit_model)
    sim = lambda th, t: base(th, t) / scale  # noqa: E731
    best = np.inf
    for k in KSTARTS:
        mle = fit_least_squares(T, data, channels, 1.0, np.log(RATES + [k]), simulate=sim, max_nfev=500)
        best = min(best, float(np.sqrt(np.mean((observed(sim(mle, T), channels) - data) ** 2))))
    return 100.0 * best                              # percent of peak (data is normalized to peak ~1)


def main():
    log("STRUCTURAL DISTINGUISHABILITY of M2 (transcriptional) vs M3 (translational) autoregulation.")
    log("Each entry = how well the WRONG mechanism mimics the right one's noise-free output (RMS residual,")
    log("percent of peak). Below the measurement noise => indistinguishable; above => distinguishable.")
    log("Deterministic, so stable. Compare against typical noise of ~1-10%.")
    log("=" * 84)
    grid = {}
    for ch, lbl in [((1,), "protein-only"), ((0, 1), "mRNA+protein")]:
        log(f"\n[{lbl}]")
        log(f"  {'1/kappa':>8} | {'M2 mimics M3':>14} | {'M3 mimics M2':>14} | {'min (best mimic)':>16}")
        for kappa in KAPPAS:
            c2 = M.make_dimensional_simulator("M2")(np.log(RATES + [kappa]), T)
            c3 = M.make_dimensional_simulator("M3")(np.log(RATES + [kappa]), T)
            r_23 = mimic_residual("M2", c3, ch)       # transcriptional model imitating translational data
            r_32 = mimic_residual("M3", c2, ch)
            grid[(lbl, kappa)] = (r_23, r_32)
            log(f"  {INV[kappa]:>8} | {r_23:>13.3f}% | {r_32:>13.3f}% | {min(r_23, r_32):>15.3f}%")

    log("\n" + "=" * 84)
    log("Reading: the residual is how distinguishable the mechanisms are. If it stays below realistic")
    log("noise (~1-10%) the mechanisms cannot be told apart from that channel; if it rises above, they can.")

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    inv = [INV[k] for k in KAPPAS]
    for ch_lbl, style in [("protein-only", "o-"), ("mRNA+protein", "s-")]:
        # the easier-to-mimic direction (min) = the harder case for telling them apart
        y = [min(grid[(ch_lbl, k)]) for k in KAPPAS]
        ax.plot(inv, y, style, label=ch_lbl)
    for noise, c in [(1, "0.7"), (3, "0.5"), (10, "0.3")]:
        ax.axhline(noise, color=c, lw=0.8, ls=":")
        ax.text(inv[0], noise * 1.05, f"{noise}% noise", fontsize=7, color=c)
    ax.set_yscale("log")
    ax.set_xlabel("feedback strength 1/kappa (molecules; left = stronger)")
    ax.set_ylabel("best wrong-mechanism mimic residual (% of peak, log)")
    ax.set_title("Distinguishability of transcriptional vs translational autoregulation")
    ax.legend(fontsize=8)
    fig.tight_layout()
    f = os.path.join(RESULTS, "fig13_structural_distinguishability.png")
    fig.savefig(f, dpi=140)
    log(f"wrote {os.path.relpath(f)}")
    out = os.path.join(RESULTS, "structural_distinguishability_summary.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    log(f"wrote {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
