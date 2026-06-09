"""Phase 3 (real data): load and inspect the Frohlich et al. 2018 single-cell GFP translation data.

This is real single-cell fluorescence data for an mRNA->protein system (chemically modified mRNA encoding
eGFP, transfected as a bolus at t=0; GFP measured by time-lapse, mRNA latent). It is the dataset used in
the published structural-identifiability study Pieschner, Hasenauer & Fuchs (2022), so it lets us fit a
small mRNA->protein ODE to real data and check the identifiability result against a known answer.

Source: Zenodo 10.5281/zenodo.1228899 (code.zip), file code/project/data/20160427_mean_eGFP.xlsx.
Licensed CC-BY-SA-4.0 (Frohlich et al. 2018, npj Syst Biol Appl, DOI 10.1038/s41540-018-0079-7). The
raw zip lives in data/raw (gitignored); this script downloads it if missing.

Run: python experiments/05_real_data/load_and_inspect.py
"""
from __future__ import annotations

import os
import zipfile

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(__file__)
RAW = os.path.join(HERE, "..", "..", "data", "raw")
RESULTS = os.path.join(HERE, "..", "..", "results")
os.makedirs(RESULTS, exist_ok=True)
ZIP = os.path.join(RAW, "frohlich_code.zip")
URL = "https://zenodo.org/records/1228899/files/code.zip?download=1"
MEMBER = "code/project/data/20160427_mean_eGFP.xlsx"


def load_traces():
    """Return (t_seconds, traces) where traces has shape (n_timepoints, n_cells). Downloads if needed."""
    out = os.path.join(RAW, "frohlich", MEMBER)
    if not os.path.exists(out):
        if not os.path.exists(ZIP):
            import urllib.request
            print("downloading Frohlich data from Zenodo (~100 MB)...")
            urllib.request.urlretrieve(URL, ZIP)
        with zipfile.ZipFile(ZIP) as z:
            z.extract(MEMBER, os.path.join(RAW, "frohlich"))
    arr = pd.read_excel(out, header=None).values.astype(float)
    return arr[:, 0], arr[:, 1:]


def main():
    t, traces = load_traces()
    mean = np.nanmean(traces, axis=1)
    sd = np.nanstd(traces, axis=1)
    th = t / 3600.0
    print(f"timepoints = {len(t)}, cells = {traces.shape[1]}, "
          f"t = {t[0]:.0f}..{t[-1]:.0f} s ({th[-1]:.1f} h), dt = {t[1] - t[0]:.0f} s")
    print(f"mean signal: min {mean.min():.2f} at {th[np.argmin(mean)]:.1f} h, "
          f"max {mean.max():.2f} at {th[np.argmax(mean)]:.1f} h, end {mean[-1]:.2f}")
    print("shape is the translation signature if it rises to a peak then decays (mRNA bolus depletes).")

    fig, ax = plt.subplots(figsize=(7, 4))
    step = max(1, traces.shape[1] // 12)
    for i in range(0, traces.shape[1], step):
        ax.plot(th, traces[:, i], color="0.8", lw=0.5)
    ax.plot(th, mean, "b", lw=2, label=f"mean of {traces.shape[1]} cells")
    ax.fill_between(th, mean - sd, mean + sd, color="b", alpha=0.15, label="+/- 1 sd")
    ax.set_xlabel("time (h)")
    ax.set_ylabel("GFP fluorescence (a.u.)")
    ax.set_title("Frohlich et al. 2018 single-cell eGFP translation kinetics (real data)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    f = os.path.join(RESULTS, "fig07_frohlich_data.png")
    fig.savefig(f, dpi=140)
    print(f"wrote {os.path.relpath(f)}")


if __name__ == "__main__":
    main()
