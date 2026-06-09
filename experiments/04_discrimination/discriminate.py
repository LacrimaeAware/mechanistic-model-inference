"""Phase 4: model discrimination. Can the models be told apart from data, and under what measurement?

Method: generate noisy data from a known "truth" model, fit each candidate model by maximum likelihood,
and compare with AIC and BIC (Gaussian likelihood, known sigma; the additive constant is identical
across candidates on the same data, so only 2k + RSS/sigma^2 and k ln N + RSS/sigma^2 are compared).
Lower is better; the reported delta is candidate minus best.

Two scenarios, each testing a Phase 1 hypothesis:
  A. truth = M3, candidates {M1, M3}. Hypothesis: the mRNA channel cannot separate M1 from M3 (their
     mRNA dynamics are identical), so only protein observation distinguishes them.
  B. truth in {M2, M3} at symmetric feedback, candidates {M2, M3}. Hypothesis: the mRNA channel
     separates them (only M2 regulates mRNA), while protein alone is weaker because their protein
     steady states coincide under symmetric parameters.

Run: python experiments/04_discrimination/discriminate.py
"""
from __future__ import annotations

import os
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from mechanistic_inference import models as M  # noqa: E402
from mechanistic_inference import fit_least_squares, observed  # noqa: E402

RESULTS = os.path.join(os.path.dirname(__file__), "..", "..", "results")
os.makedirs(RESULTS, exist_ok=True)
T = np.linspace(0.0, 300.0, 80)
SIGMA = 0.03
lines = []

# true dimensional parameters per model (k_m, d_m, k_p, d_p[, kappa]); kappa = 0.02 (1/kappa = 50)
TRUTH = {
    "M1": np.log([1.04, 0.35, 7.70, 0.02]),
    "M2": np.log([1.04, 0.35, 7.70, 0.02, 0.02]),
    "M3": np.log([1.04, 0.35, 7.70, 0.02, 0.02]),
}
START = {  # generic starting guess for fitting a candidate (rates from Table 1, kappa default)
    "M1": np.log([1.04, 0.35, 7.70, 0.02]),
    "M2": np.log([1.04, 0.35, 7.70, 0.02, 0.02]),
    "M3": np.log([1.04, 0.35, 7.70, 0.02, 0.02]),
}


def record(s):
    print(s); lines.append(s)


def fit_aic_bic(candidate, data, ch, scale):
    base = M.make_dimensional_simulator(candidate)
    sim_norm = lambda th, t: base(th, t) / scale  # noqa: E731
    mle = fit_least_squares(T, data, ch, SIGMA, START[candidate], simulate=sim_norm, max_nfev=300)
    pred = observed(sim_norm(mle, T), ch)
    rss = float(np.sum((pred - np.asarray(data)) ** 2))
    k, N = len(START[candidate]), len(data)
    return {"aic": 2 * k + rss / SIGMA ** 2, "bic": k * np.log(N) + rss / SIGMA ** 2, "rss": rss}


def scenario(truth, candidates, schemes, seed):
    base = M.make_dimensional_simulator(truth)
    clean = base(TRUTH[truth], T)
    scale = np.array([clean[:, 0].max(), clean[:, 1].max()])
    record(f"truth = {truth};  candidates = {candidates}")
    for ch, lbl in schemes:
        rng = np.random.default_rng(seed)
        data = (clean / scale)[:, list(ch)].ravel() + rng.normal(scale=SIGMA, size=len(ch) * T.size)
        res = {c: fit_aic_bic(c, data, ch, scale) for c in candidates}
        best = min(candidates, key=lambda c: res[c]["aic"])
        parts = []
        for c in candidates:
            dA = res[c]["aic"] - res[best]["aic"]
            parts.append(f"{c} dAIC={dA:7.1f}")
        record(f"  [{lbl:12s}] winner={best}  " + "  ".join(parts))
    record("")


def main():
    record("SCENARIO A: truth M3, can M1 vs M3 be told apart, and from which channel?")
    record("=" * 70)
    scenario("M3", ["M1", "M3"],
             [((0,), "mRNA-only"), ((1,), "protein-only"), ((0, 1), "mRNA+protein")], seed=1)
    record("Expectation: mRNA-only cannot separate (identical mRNA, so M1 wins on parsimony);")
    record("protein observation reveals M3. dAIC near 0 (favoring the simpler M1) under mRNA-only;")
    record("large dAIC favoring M3 under protein.")
    record("")

    record("SCENARIO B: truth M2 or M3 (symmetric feedback), can M2 vs M3 be told apart?")
    record("=" * 70)
    scenario("M2", ["M2", "M3"],
             [((0,), "mRNA-only"), ((1,), "protein-only"), ((0, 1), "mRNA+protein")], seed=2)
    scenario("M3", ["M2", "M3"],
             [((0,), "mRNA-only"), ((1,), "protein-only"), ((0, 1), "mRNA+protein")], seed=3)
    record("Expectation: the mRNA channel separates M2 from M3 (only M2 regulates mRNA); protein alone")
    record("is weaker because the protein steady states coincide under symmetric parameters.")
    record("")

    out = os.path.join(RESULTS, "discrimination_summary.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    record(f"wrote {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
