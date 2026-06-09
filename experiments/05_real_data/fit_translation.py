"""Phase 3 (real data): fit the translation model to the Frohlich mean GFP trace and check identifiability.

PROVISIONAL. This is a first pass and has not been independently re-verified. It (a) fits the closed-form
translation model to the real mean GFP curve with multi-start, (b) checks the fit visually and by
residual, and (c) asks whether our Fisher-rank tool finds the published structural non-identifiability:
the amplitude is scale*k*m0, so from GFP alone those three factors are not separable (two null
directions), which the Pieschner et al. 2022 analysis establishes for this model.

Run: python experiments/05_real_data/fit_translation.py
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
from mechanistic_inference import fit_least_squares, fisher_information, observed  # noqa: E402
from load_and_inspect import load_traces  # noqa: E402

RESULTS = os.path.join(os.path.dirname(__file__), "..", "..", "results")
os.makedirs(RESULTS, exist_ok=True)
lines = []


def log(s):
    print(s, flush=True); lines.append(s)


def main():
    t_s, traces = load_traces()
    th = t_s / 3600.0                       # hours
    y = np.nanmean(traces, axis=1)          # mean GFP trace
    log(f"mean trace: {len(th)} points, 0..{th[-1]:.1f} h, signal {y.min():.0f}..{y.max():.0f} a.u.")

    # multi-start fit of simulate_translation: theta = log[amplitude, delta, gamma, t0, offset]
    starts = []
    for delta in (0.05, 0.1, 0.2):
        for gamma in (0.02, 0.08, 0.3):
            if abs(delta - gamma) < 1e-6:
                continue
            starts.append(np.log([2.0e3, delta, gamma, 0.2, max(y[0], 1.0)]))
    best_rss, best = np.inf, None
    for s in starts:
        mle = fit_least_squares(th, y, (0,), 1.0, s, simulate=M.simulate_translation,
                                log_bound=8.0, max_nfev=600)
        rss = float(np.sum((observed(M.simulate_translation(mle, th), (0,)) - y) ** 2))
        if rss < best_rss:
            best_rss, best = rss, mle
    amp, delta, gamma, t0, offset = np.exp(best)
    resid_sd = np.sqrt(best_rss / len(y))
    log(f"best fit: amplitude={amp:.1f}, delta={delta:.3f}/h, gamma={gamma:.3f}/h, t0={t0:.2f} h, "
        f"offset={offset:.1f}; residual sd={resid_sd:.1f} a.u. ({100 * resid_sd / y.max():.1f}% of peak)")
    log("(delta and gamma are the two decay rates; only the unordered pair is identifiable, the labels"
        " are not, which is a known symmetry of this model.)")

    fit_curve = observed(M.simulate_translation(best, th), (0,))
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(th, y, "k.", ms=3, label="real mean GFP (1079 cells)")
    ax.plot(th, fit_curve, "r-", lw=2, label="translation-model fit")
    ax.set_xlabel("time (h)"); ax.set_ylabel("GFP fluorescence (a.u.)")
    ax.set_title("Frohlich mean trace: translation-model fit (provisional)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    f = os.path.join(RESULTS, "fig09_translation_fit.png")
    fig.savefig(f, dpi=140)
    log(f"wrote {os.path.relpath(f)}")

    # Identifiability check: Fisher rank in the SEPARATED parameterization [scale, k, m0, delta, gamma, t0, offset]
    cube = amp ** (1.0 / 3.0)               # split the amplitude so scale*k*m0 = amp
    theta_full = np.log([cube, cube, cube, delta, gamma, t0, offset])
    names = ["scale", "k", "m0", "delta", "gamma", "t0", "offset"]
    F = fisher_information(theta_full, th, (0,), max(resid_sd, 1e-6), simulate=M.simulate_translation_full)
    w, V = np.linalg.eigh(F)
    ev = np.clip(w, 0, None)
    rank = int(np.sum(ev > 1e-9 * (ev.max() + 1e-30)))
    log("")
    log(f"Fisher rank (separated 7-parameter model): {rank} of 7  "
        f"(expected 5: scale, k, m0 enter only as a product, so 2 null directions)")
    for j in range(min(2, V.shape[1])):           # the two smallest-eigenvalue directions
        d = ", ".join(f"{n}:{v:+.2f}" for n, v in zip(names, V[:, j]) if abs(v) > 0.15)
        log(f"  null direction {j + 1} (eig {ev[j]:.1e} rel {ev[j] / ev.max():.1e}): {d}")
    log("Reading: if the two null directions live in the scale/k/m0 block, the pipeline has found the")
    log("published product non-identifiability on real data. PROVISIONAL until independently re-checked.")

    out = os.path.join(RESULTS, "translation_fit_summary.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    log(f"wrote {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
