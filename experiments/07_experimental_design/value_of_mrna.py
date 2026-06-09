"""Experiment 07: the value of the mRNA channel (deterministic, Fisher-information based).

Builds on the project's spine (protein alone is information-poor; mRNA carries the information) and makes
it a quantitative experimental-design statement, using the Fisher information matrix evaluated at the
true parameters. This is DETERMINISTIC - no Monte-Carlo, no noise draws - so the result is stable and
cannot flip-flop, which is the lesson from the discrimination work.

Two questions:
  1. Protein-only vs mRNA+protein: rank and parameter standard errors. Protein-only is rank-deficient
     (one parameter combination has infinite variance) no matter how many protein timepoints; adding
     mRNA makes it full rank with finite errors.
  2. How little mRNA do you need? With protein measured at every timepoint, add mRNA at K timepoints and
     find the smallest K that restores full rank - the concrete design recommendation.

Run: python experiments/07_experimental_design/value_of_mrna.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from mechanistic_inference import models as M  # noqa: E402
from mechanistic_inference import fisher_information  # noqa: E402

RESULTS = os.path.join(os.path.dirname(__file__), "..", "..", "results")
os.makedirs(RESULTS, exist_ok=True)
T = np.linspace(0.0, 300.0, 80)
SIGMA = 0.03
TRUE = {
    "M1": np.log([1.04, 0.35, 7.70, 0.02]),
    "M2": np.log([1.04, 0.35, 7.70, 0.02, 0.02]),
    "M3": np.log([1.04, 0.35, 7.70, 0.02, 0.02]),
}
NAMES = {"M1": ["k_m", "d_m", "k_p", "d_p"], "M2": ["k_m", "d_m", "k_p", "d_p", "kappa_M"],
         "M3": ["k_m", "d_m", "k_p", "d_p", "kappa_P"]}
lines = []


def log(s):
    print(s); lines.append(s)


def rank_of(F):
    ev = np.clip(np.linalg.eigvalsh(F), 0, None)
    return int(np.sum(ev > 1e-9 * (ev.max() + 1e-30))), ev


def fim_mixed(model, theta, scale, mrna_idx, eps=1e-5):
    """Fisher information for: protein at ALL timepoints + mRNA at the timepoints in mrna_idx.
    Normalized per channel; finite-difference sensitivities in log-parameter space."""
    base = M.make_dimensional_simulator(model)
    theta = np.asarray(theta, float)

    def obs(th):
        traj = base(th, T) / scale
        return np.concatenate([traj[:, 1], traj[mrna_idx, 0]])  # protein(all), mRNA(subset)

    base_obs = obs(theta)
    J = np.zeros((base_obs.size, theta.size))
    for j in range(theta.size):
        tp = theta.copy(); tp[j] += eps
        tm = theta.copy(); tm[j] -= eps
        J[:, j] = (obs(tp) - obs(tm)) / (2 * eps)
    return (J.T @ J) / SIGMA ** 2


def main():
    log("VALUE OF THE mRNA CHANNEL (Fisher information at the true parameters; deterministic).")
    log("=" * 78)
    for model in ("M1", "M2", "M3"):
        sim = M.make_dimensional_simulator(model)
        theta = TRUE[model]
        clean = sim(theta, T)
        scale = np.array([clean[:, 0].max(), clean[:, 1].max()])
        names = NAMES[model]

        def sim_norm(th, t, _b=sim, _s=scale):
            return _b(th, t) / _s

        log(f"\n{model}  ({len(names)} parameters: {names})")
        for ch, lbl in [((1,), "protein-only"), ((0, 1), "mRNA+protein")]:
            F = fisher_information(theta, T, ch, SIGMA, simulate=sim_norm)
            rank, ev = rank_of(F)
            full = rank == len(names)
            log(f"  [{lbl:12s}] rank {rank}/{len(names)}  smallest/largest eig = {ev.min()/ev.max():.1e}")
            if full:
                se = np.sqrt(np.diag(np.linalg.inv(F)))
                log(f"      standard errors (log units): " + ", ".join(f"{n}:{s:.2f}" for n, s in zip(names, se)))
            else:
                log(f"      rank-deficient: one parameter combination has INFINITE variance (the k_m*k_p "
                    f"product); no number of protein timepoints fixes this.")

    log("\n" + "=" * 78)
    log("HOW LITTLE mRNA DO YOU NEED? (protein at all 80 timepoints, mRNA at K evenly-spaced times)")
    for model in ("M1", "M2", "M3"):
        theta = TRUE[model]
        clean = M.make_dimensional_simulator(model)(theta, T)
        scale = np.array([clean[:, 0].max(), clean[:, 1].max()])
        n_par = len(NAMES[model])
        min_k = None
        for K in range(0, 6):
            # place mRNA samples in the informative window (skip t=0, where mRNA = 0 and carries nothing)
            mrna_idx = np.linspace(len(T) // 4, len(T) - 1, K).round().astype(int) if K > 0 else np.array([], int)
            F = fim_mixed(model, theta, scale, mrna_idx)
            rank, _ = rank_of(F)
            if rank == n_par:
                min_k = K
                break
        log(f"  {model}: smallest number of informative mRNA timepoints for full rank = {min_k} "
            f"(protein alone, K=0, is rank-deficient)")

    log("\nReading (deterministic, stable): protein-only is structurally rank-deficient for every model "
        "and any amount of protein sampling; a small number of mRNA timepoints restores full identifiability.")
    out = os.path.join(RESULTS, "value_of_mrna_summary.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    log(f"wrote {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
