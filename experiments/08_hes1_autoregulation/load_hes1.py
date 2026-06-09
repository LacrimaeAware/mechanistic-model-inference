"""Experiment 08: load and characterize real Hes1 single-cell oscillation data before modeling.

Hes1 is a real negative-autoregulation gene (it represses its own transcription), in the family of the
source paper's transcriptional model. Data: GPosc repo (Phillips et al. 2017), Hes1_example.csv -
12 single-cell fluorescent-reporter traces, 0.5 h sampling over ~70 h, with background columns.

Each cell is tracked for a different length of time (its column goes NaN when the cell is lost), so every
cell is handled on its own valid (non-NaN) window. This script only looks at the data: background-
subtract, plot, and estimate whether and at what period the traces oscillate - the dynamics decide what
model is appropriate (a plain non-delay model only damps; sustained oscillation needs a delay).

Run: python experiments/08_hes1_autoregulation/load_hes1.py
"""
from __future__ import annotations

import os
import warnings

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
HERE = os.path.dirname(__file__)
RAW = os.path.join(HERE, "..", "..", "data", "raw", "hes1", "Hes1_example.csv")
RESULTS = os.path.join(HERE, "..", "..", "results")
os.makedirs(RESULTS, exist_ok=True)
lines = []


def log(s):
    print(s); lines.append(s)


def cell_period(tv, yv):
    """Dominant oscillation period (h) of a detrended single-cell trace, or None."""
    if len(tv) < 20:
        return None
    osc = yv - np.polyval(np.polyfit(tv, yv, 3), tv)         # remove slow trend
    ac = np.correlate(osc - osc.mean(), osc - osc.mean(), mode="full")[len(osc) - 1:]
    ac = ac / ac[0]
    dt = tv[1] - tv[0]
    for k in range(2, len(ac) - 1):
        if ac[k] > ac[k - 1] and ac[k] > ac[k + 1] and ac[k] > 0.1:
            return k * dt
    return None


def main():
    df = pd.read_csv(RAW)
    t = df["Time (h)"].values
    cell_cols = [c for c in df.columns if c.startswith("Cell")]
    bg_cols = [c for c in df.columns if c.startswith("Background")]
    cells = df[cell_cols].values
    bg = np.nanmean(df[bg_cols].values, axis=1)              # per-timepoint background, ignoring NaN
    log(f"{len(t)} timepoints, {t[0]}..{t[-1]} h, dt={t[1]-t[0]} h; {len(cell_cols)} cells.")

    traces, lengths, periods = [], [], []
    for j in range(cells.shape[1]):
        m = ~np.isnan(cells[:, j]) & ~np.isnan(bg)
        tv, yv = t[m], cells[m, j] - bg[m]                  # background-subtracted, valid window only
        traces.append((tv, yv))
        if len(tv) > 5:
            lengths.append(tv[-1] - tv[0])
        p = cell_period(tv, yv)
        if p:
            periods.append(p)
    log(f"tracked durations: median {np.median(lengths):.1f} h [{np.min(lengths):.1f}-{np.max(lengths):.1f}].")
    if periods:
        log(f"dominant period (autocorrelation of detrended traces): median {np.median(periods):.1f} h "
            f"[{np.percentile(periods,25):.1f}-{np.percentile(periods,75):.1f}], detected in "
            f"{len(periods)}/{cells.shape[1]} cells. Hes1 is reported to oscillate ~2-3 h.")
    else:
        log("no clear oscillation period detected by autocorrelation.")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for tv, yv in traces:
        axes[0].plot(tv, yv, lw=0.8, alpha=0.75)
    axes[0].set_xlabel("time (h)"); axes[0].set_ylabel("Hes1 reporter (bg-subtracted)")
    axes[0].set_title(f"{len(traces)} real Hes1 single cells (own valid windows)")
    off = 0
    for tv, yv in traces[:4]:
        osc = yv - np.polyval(np.polyfit(tv, yv, 3), tv)
        axes[1].plot(tv, osc + off, lw=0.9); off += 400
    axes[1].set_xlabel("time (h)"); axes[1].set_ylabel("detrended (offset per cell)")
    axes[1].set_title("oscillation component (detrended)")
    fig.tight_layout()
    f = os.path.join(RESULTS, "fig12_hes1_data.png")
    fig.savefig(f, dpi=140)
    log(f"wrote {os.path.relpath(f)}")
    out = os.path.join(RESULTS, "hes1_data_summary.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    log(f"wrote {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
