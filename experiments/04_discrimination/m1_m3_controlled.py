"""Controlled, definitive check of the contested cell: M1 vs M3 discrimination from protein only.

Motivation: the single-run discrimination and the replication gave different impressions of this cell,
partly because they used different grids and different optimizers. This script holds EVERYTHING fixed
(grid, optimizer, noise, feedback strength) and reports the full distribution of delta-AIC over many
noise draws, so the answer is a distribution, not a flip-flopping point estimate.

Built-in correctness check: M1 is M3 with kappa_P = 0, so M1 is nested in M3 and the more flexible M3
cannot fit worse. Therefore delta-AIC = AIC(M1) - AIC(M3) = -2 + (RSS_M1 - RSS_M3)/sigma^2 must be
>= -2 for every draw. Any draw below -2 is an optimizer failure (the M3 fit did not reach the M1 fit),
not a real result; the script counts those.

Run: python -u experiments/04_discrimination/m1_m3_controlled.py
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

# fixed settings (stated once, used for every fit)
T = np.linspace(0.0, 300.0, 80)
SIGMA = 0.03
N = 60
MAX_NFEV = 500              # generous, so the nested M3 fit reaches the M1 fit
KAPPA = 0.02               # 1/kappa = 50 molecules (feedback strength), fixed
CH = (1,)                  # protein only
TRUE_M3 = np.log([1.04, 0.35, 7.70, 0.02, KAPPA])
START = {"M1": np.log([1.04, 0.35, 7.70, 0.02]), "M3": TRUE_M3}
lines = []


def log(s):
    print(s, flush=True)
    lines.append(s)


def best_fit(model, data, scale, starts):
    """Multi-start fit: try every start, return the (lowest RSS, its theta). Multi-start is needed
    because the kappa direction is flat near zero, so a single start gets stuck (the bug this check
    found)."""
    base = M.make_dimensional_simulator(model)
    sim = lambda th, t: base(th, t) / scale  # noqa: E731
    best_rss, best_mle = np.inf, None
    for s in starts:
        mle = fit_least_squares(T, data, CH, SIGMA, s, simulate=sim, max_nfev=MAX_NFEV)
        rss = float(np.sum((observed(sim(mle, T), CH) - np.asarray(data)) ** 2))
        if rss < best_rss:
            best_rss, best_mle = rss, mle
    return best_rss, best_mle


def main():
    log(f"M1 vs M3, protein-only, FIXED settings: grid={len(T)}, sigma={SIGMA}, 1/kappa=50, "
        f"least-squares max_nfev={MAX_NFEV}, N={N} draws.")
    log("delta-AIC = AIC(M1) - AIC(M3); positive favors the true model M3. Floor is -2 (nesting).")
    log("=" * 80)
    base3 = M.make_dimensional_simulator("M3")
    clean = base3(TRUE_M3, T)
    scale = np.array([clean[:, 0].max(), clean[:, 1].max()])

    deltas, below_floor = [], 0
    examples = []
    for d in range(N):
        rng = np.random.default_rng(d)
        data = (clean / scale)[:, list(CH)].ravel() + rng.normal(scale=SIGMA, size=T.size)
        rss1, mle1 = best_fit("M1", data, scale, [START["M1"]])
        # M3 multi-start: the M1 fit + tiny kappa (guarantees the nesting floor) AND several
        # kappa-active starts (so a real feedback signature can be found despite the flat-near-zero gradient)
        m3_starts = [np.append(mle1, np.log(k)) for k in (1e-3, 5e-3, 2e-2, 8e-2, 3e-1)]
        rss3, _ = best_fit("M3", data, scale, m3_starts)
        aic1, aic3 = 2 * 4 + rss1 / SIGMA ** 2, 2 * 5 + rss3 / SIGMA ** 2
        delta = aic1 - aic3
        deltas.append(delta)
        if delta < -2.0 - 1e-6:
            below_floor += 1
        if d < 5:
            examples.append((d, rss1, rss3, delta))
    deltas = np.array(deltas)

    log("five individual draws (to show the spread a single run hides):")
    for d, rss1, rss3, delta in examples:
        log(f"  seed {d}: RSS_M1={rss1:.4f}  RSS_M3={rss3:.4f}  delta-AIC={delta:+.2f}  "
            f"-> {'M3' if delta > 0 else 'M1'} wins")
    log("")
    q = np.percentile(deltas, [5, 25, 50, 75, 95])
    log(f"delta-AIC distribution over {N} draws:")
    log(f"  mean {deltas.mean():+.2f}, std {deltas.std():.2f}, min {deltas.min():+.2f}, max {deltas.max():+.2f}")
    log(f"  quantiles 5/25/50/75/95 = {q[0]:+.1f} / {q[1]:+.1f} / {q[2]:+.1f} / {q[3]:+.1f} / {q[4]:+.1f}")
    rate = float(np.mean(deltas > 0))
    near = float(np.mean((deltas > 0) & (deltas < 2)))
    log(f"  M3 (true) selected: {rate:.0%}   (of which within 2 AIC of M1, i.e. weak: {near:.0%})")
    log(f"  correctness check: draws below the -2 nesting floor = {below_floor}/{N} "
        f"({'OK' if below_floor == 0 else 'OPTIMIZER FAILURES PRESENT'})")
    log("")
    log("Conclusion: a single draw is one sample from this whole distribution. The earlier single-run")
    log("number (delta ~ +0.6) and the replication mean are both samples/summaries of THIS spread; the")
    log("honest statement is the distribution and the selection rate, not any one delta-AIC.")

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(deltas, bins=20, color="steelblue", alpha=0.8)
    ax.axvline(0, color="k", lw=1, ls="-", label="decision boundary")
    ax.axvline(-2, color="r", lw=1, ls=":", label="nesting floor (-2)")
    ax.axvline(deltas.mean(), color="orange", lw=1.5, ls="--", label=f"mean {deltas.mean():+.1f}")
    ax.set_xlabel("delta-AIC = AIC(M1) - AIC(M3)  (positive favors true M3)")
    ax.set_ylabel(f"count of {N} noise draws")
    ax.set_title(f"M1 vs M3 protein-only: true model selected {rate:.0%} of draws")
    ax.legend(fontsize=8)
    fig.tight_layout()
    f = os.path.join(RESULTS, "fig08_m1_m3_controlled.png")
    fig.savefig(f, dpi=140)
    log(f"wrote {os.path.relpath(f)}")
    out = os.path.join(RESULTS, "m1_m3_controlled_summary.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    log(f"wrote {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
