"""Phase 3 (real data), the tangible version: recover translation-kinetics parameters for real single cells.

Fit the closed-form translation model to individual Frohlich single-cell GFP traces (a sample of them)
and report the recovered, identifiable parameters and their cell-to-cell spread. This is the concrete
"fit real data -> here are the parameters" deliverable. The model is analytic, so fits are fast.

Only identifiable quantities are reported: the amplitude (scale*k*m0 product), the two decay rates as an
UNORDERED pair (fast and slow, since their labels are not identifiable), the onset time, and the offset.

Run: python experiments/05_real_data/fit_single_cells.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))
from mechanistic_inference import models as M  # noqa: E402
from mechanistic_inference import fit_least_squares, observed  # noqa: E402
from load_and_inspect import load_traces  # noqa: E402

RESULTS = os.path.join(os.path.dirname(__file__), "..", "..", "results")
os.makedirs(RESULTS, exist_ok=True)
lines = []


def log(s):
    print(s, flush=True); lines.append(s)


def fit_cell(th_hours, y):
    """Multi-start fit of the translation model to one cell. Returns (params, rms_pct) or None."""
    if y.max() - y.min() < 50:          # skip cells with essentially no signal
        return None
    best = None
    for delta in (0.05, 0.1, 0.2):
        for gamma in (0.02, 0.08, 0.3):
            if abs(delta - gamma) < 1e-6:
                continue
            s = np.log([max(y.max(), 1.0), delta, gamma, 0.2, max(y[0], 1.0)])
            mle = fit_least_squares(th_hours, y, (0,), 1.0, s, simulate=M.simulate_translation,
                                    log_bound=8.0, max_nfev=400)
            rms = float(np.sqrt(np.mean((observed(M.simulate_translation(mle, th_hours), (0,)) - y) ** 2)))
            if best is None or rms < best[1]:
                best = (np.exp(mle), rms)
    amp, delta, gamma, t0, offset = best[0]
    return np.array([amp, max(delta, gamma), min(delta, gamma), t0, offset]), 100 * best[1] / y.max()


def main():
    t_s, traces = load_traces()
    th = t_s / 3600.0
    rng = np.random.default_rng(0)
    idx = rng.choice(traces.shape[1], size=200, replace=False)
    log(f"Fitting the translation model to {len(idx)} randomly sampled real single cells (of {traces.shape[1]}).")

    rows, kept, bad = [], 0, 0
    for i in idx:
        r = fit_cell(th, traces[:, i])
        if r is None:
            continue
        params, rms = r
        if rms < 8.0:                    # keep reasonable fits (RMS < 8% of that cell's peak)
            rows.append(params); kept += 1
        else:
            bad += 1
    A = np.array(rows)
    log(f"kept {kept} good fits (RMS < 8% of peak); discarded {bad} poor fits and "
        f"{len(idx) - kept - bad} low-signal cells.")
    log("")
    names = ["amplitude", "fast rate (1/h)", "slow rate (1/h)", "onset t0 (h)", "offset (a.u.)"]
    log("Recovered parameters across cells (median [25th-75th percentile]):")
    for j, n in enumerate(names):
        q = np.percentile(A[:, j], [25, 50, 75])
        log(f"  {n:18s}: {q[1]:8.3g}  [{q[0]:.3g} - {q[2]:.3g}]")
    log("")
    log("Half-lives from the rates: fast = ln2/rate, slow = ln2/rate (hours):")
    log(f"  fast-rate half-life median {np.log(2)/np.median(A[:,1]):.2f} h; "
        f"slow-rate half-life median {np.log(2)/np.median(A[:,2]):.2f} h")

    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    for ax, j, n in zip(axes.ravel(), [0, 1, 2, 3], names):
        vals = A[:, j]
        ax.hist(vals, bins=25, color="steelblue", alpha=0.85)
        ax.axvline(np.median(vals), color="orange", lw=1.5, ls="--", label=f"median {np.median(vals):.3g}")
        ax.set_xlabel(n); ax.set_ylabel("cells"); ax.legend(fontsize=8)
    fig.suptitle(f"Recovered translation-kinetics parameters for {kept} real single cells (Frohlich 2018)")
    fig.tight_layout()
    f = os.path.join(RESULTS, "fig11_single_cell_params.png")
    fig.savefig(f, dpi=140)
    log(f"wrote {os.path.relpath(f)}")
    out = os.path.join(RESULTS, "single_cell_params_summary.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    log(f"wrote {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
