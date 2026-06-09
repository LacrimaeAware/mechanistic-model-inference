"""Reproduce the stochastic noise ordering of the three negative-autoregulation models.

Tests the paper's claim that, by the protein coefficient of variation, transcriptional autoregulation
is noisiest, the simplistic model next, and translational autoregulation least. Uses the exact SSA in
mechanistic_inference.stochastic. See that module for the propensity caveat (the Hill-in-propensity
formulation is a modeling choice, so the qualitative ordering is the target, not the exact CV values).

Run: python experiments/01_reproduce_autoregulation_models/reproduce_stochastic.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from mechanistic_inference import stochastic as S  # noqa: E402

RESULTS = os.path.join(os.path.dirname(__file__), "..", "..", "results")
os.makedirs(RESULTS, exist_ok=True)


def collect(seed, **params):
    """Pool protein/mRNA counts over an ensemble at the sample times, for CV and histograms."""
    rng = np.random.default_rng(seed)
    st = [250, 300, 350, 400, 450, 500, 550, 600]
    Ms, Ps = [], []
    for _ in range(params.pop("n_runs")):
        s = S.gillespie_nar(sample_times=st, rng=rng, **params)
        Ms.append(s[:, 0]); Ps.append(s[:, 1])
    return np.concatenate(Ms), np.concatenate(Ps)


def main():
    lines = []

    def record(s):
        print(s); lines.append(s)

    base = dict(alpha_M=1.04, beta_M=0.35, beta_P=0.02, n_M=2.0, n_P=2.0, S_a=0.0, M0=3, P0=400, n_runs=300)
    kap = 1.0 / 300.0  # 1/kappa = 300 molecules

    runs = {
        "M1 simplistic": dict(seed=1, alpha_P=2.94, **base),                  # alpha_P reset so the mean protein is comparable
        "M2 transcriptional": dict(seed=2, alpha_P=7.70, kappa_M=kap, **base),
        "M3 translational": dict(seed=3, alpha_P=7.70, kappa_P=kap, **base),
    }

    stats = {}
    dists = {}
    for name, params in runs.items():
        M, P = collect(**params)
        cv = lambda x: x.std() / x.mean()
        stats[name] = (P.mean(), cv(P), M.mean(), cv(M))
        dists[name] = P
        record(f"{name:20s} mean protein={P.mean():6.1f}  CV_protein={cv(P):.3f}  "
               f"mean mRNA={M.mean():5.2f}  CV_mRNA={cv(M):.3f}")

    cvP = {k: v[1] for k, v in stats.items()}
    record(f"protein-CV ordering M2 > M1 > M3: "
           f"{cvP['M2 transcriptional'] > cvP['M1 simplistic'] > cvP['M3 translational']}")
    record(f"  (M2={cvP['M2 transcriptional']:.3f} > M1={cvP['M1 simplistic']:.3f} > "
           f"M3={cvP['M3 translational']:.3f})")
    record("paper reported CV_protein: M1 0.14, M2 0.16, M3 0.10; CV_mRNA: M1 0.56, M2 0.89, M3 0.56")

    fig, ax = plt.subplots(figsize=(6, 4))
    for name, P in dists.items():
        ax.hist(P, bins=40, density=True, histtype="step", linewidth=1.6,
                label=f"{name} (CV={stats[name][1]:.3f})")
    ax.set_xlabel("protein count")
    ax.set_ylabel("probability density")
    ax.set_title("Steady-state protein distributions (1/kappa = 300 molecules, S_a = 0)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    f3 = os.path.join(RESULTS, "fig03_stochastic_cv.png")
    fig.savefig(f3, dpi=140)
    record(f"wrote {os.path.relpath(f3)}")

    out = os.path.join(RESULTS, "reproduce_01_stochastic_summary.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    record(f"wrote {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
