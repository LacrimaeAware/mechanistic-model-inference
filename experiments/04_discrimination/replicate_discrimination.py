"""Phase 4 backward verification: replicate the borderline AIC discrimination cells over noise draws.

A single delta-AIC near 0 is not evidence of "indistinguishable"; it could be one noise draw. This
repeats the cells that were borderline in the single run (discriminate.py) over N independent noise
realizations and reports the model-selection rate: how often the data-generating model is the AIC
winner. ~50% means genuinely indistinguishable; ~100% means reliably distinguishable.

Settings are deliberately cheap (coarse grid, low optimizer cap, few draws): the rate only needs to
separate "~50%" from "~100%", not high precision. Progress is flushed after each cell and every few
draws so the run can be monitored live (read the output file, or run with python -u).

Run: python -u experiments/04_discrimination/replicate_discrimination.py
"""
from __future__ import annotations

import os
import sys
import time
import warnings

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from mechanistic_inference import models as M  # noqa: E402
from mechanistic_inference import fit_mle, observed  # noqa: E402

RESULTS = os.path.join(os.path.dirname(__file__), "..", "..", "results")
os.makedirs(RESULTS, exist_ok=True)
T = np.linspace(0.0, 300.0, 30)   # coarse grid: cheap fits
SIGMA = 0.03
N = 10                            # draws per cell (enough to tell ~50% from ~100%)
MAXITER = 150
lines = []

TRUTH = {
    "M1": np.log([1.04, 0.35, 7.70, 0.02]),
    "M2": np.log([1.04, 0.35, 7.70, 0.02, 0.02]),
    "M3": np.log([1.04, 0.35, 7.70, 0.02, 0.02]),
}

# the cells that were borderline (delta-AIC < 3) in the single run
TARGETS = [
    ("M3", ["M1", "M3"], (0,), "A: M1 vs M3 (truth M3), mRNA-only"),
    ("M3", ["M1", "M3"], (1,), "A: M1 vs M3 (truth M3), protein-only"),
    ("M2", ["M2", "M3"], (1,), "B: M2 vs M3 (truth M2), protein-only"),
    ("M3", ["M2", "M3"], (0,), "B: M2 vs M3 (truth M3), mRNA-only"),
    ("M3", ["M2", "M3"], (1,), "B: M2 vs M3 (truth M3), protein-only"),
]


def log(s):
    """Record a line and flush immediately so progress is visible during the run."""
    print(s, flush=True)
    lines.append(s)


def fit_aic(candidate, data, ch, scale):
    base = M.make_dimensional_simulator(candidate)
    sim_norm = lambda th, t: base(th, t) / scale  # noqa: E731
    mle = fit_mle(T, data, ch, SIGMA, TRUTH[candidate], simulate=sim_norm, maxiter=MAXITER)
    rss = float(np.sum((observed(sim_norm(mle, T), ch) - np.asarray(data)) ** 2))
    return 2 * len(TRUTH[candidate]) + rss / SIGMA ** 2


def main():
    log(f"Replicated borderline discrimination, N={N} draws, {SIGMA:.0%} noise, grid={len(T)}, maxiter={MAXITER}.")
    log("'select rate' = fraction of draws the true model wins on AIC (~50% = indistinguishable).")
    log("=" * 78)
    t_start = time.time()
    for cid, (truth, cands, ch, label) in enumerate(TARGETS):
        log(f"[cell {cid + 1}/{len(TARGETS)}] {label}: starting...")
        t0 = time.time()
        base = M.make_dimensional_simulator(truth)
        clean = base(TRUTH[truth], T)
        scale = np.array([clean[:, 0].max(), clean[:, 1].max()])
        other = [c for c in cands if c != truth][0]
        wins, deltas = 0, []
        for d in range(N):
            rng = np.random.default_rng(cid * 1000 + d)
            data = (clean / scale)[:, list(ch)].ravel() + rng.normal(scale=SIGMA, size=len(ch) * T.size)
            a_truth, a_other = fit_aic(truth, data, ch, scale), fit_aic(other, data, ch, scale)
            wins += int(a_truth <= a_other)
            deltas.append(a_other - a_truth)
            if (d + 1) % 3 == 0 or d == N - 1:
                log(f"    {d + 1}/{N} draws ({time.time() - t0:.0f}s, running select rate {wins/(d+1):.0%})")
        rate = wins / N
        tag = "indistinguishable" if rate < 0.65 else "reliable" if rate > 0.9 else "marginal"
        log(f"[cell {cid + 1}/{len(TARGETS)}] DONE: select {truth} {rate:.0%}  "
            f"(dAIC {np.mean(deltas):+.1f} +/- {np.std(deltas):.1f}) -> {tag}  [{time.time() - t0:.0f}s]")
    log(f"all cells done in {time.time() - t_start:.0f}s")
    out = os.path.join(RESULTS, "discrimination_replication_summary.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    log(f"wrote {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
