"""Stochastic (Gillespie) simulation of the negative-autoregulation models.

This is the optional stochastic half of Phase 1: test whether the paper's reported noise ordering
(transcriptional autoregulation noisiest, simplistic next, translational least, measured by the
protein coefficient of variation) reproduces under a standard exact stochastic simulation.

The reaction set is the mass-action reading of the dimensional ODEs (Eqs. 1-2; Table-1 parameters in
molecule counts):

    R1 transcription   : M -> M+1 ,  propensity alpha_M (1 + S_a) f(P)
    R2 mRNA decay      : M -> M-1 ,  propensity beta_M M
    R3 translation     : M -> M, P+1, propensity alpha_P M g(P)
    R4 protein decay   : P -> P-1 ,  propensity beta_P P

with f(P) = 1/(1+(kappa_M P)^n_M) and g(P) = 1/(1+(kappa_P P)^n_P).

Caveat (open, not settled): the paper does not write out how the Hill repression enters the
propensity. The choice here is to evaluate the Hill factor at the current protein count and multiply
the production propensity by it, which is the common convention. A different formulation (for example
a binding-equilibrium occupancy or a Hill on a promoter-state variable) could shift the quantitative
CV values. The qualitative ordering is the testable target; the exact CVs are not claimed to match.
"""
from __future__ import annotations

import numpy as np


def gillespie_nar(*, alpha_M, alpha_P, beta_M, beta_P, kappa_M=0.0, kappa_P=0.0, n_M=2.0, n_P=2.0,
                  S_a=0.0, M0, P0, sample_times, rng):
    """One exact SSA trajectory. Returns counts (len(sample_times), 2) = columns [M, P] sampled at
    sample_times (the state is piecewise constant between reaction events)."""
    st = np.sort(np.asarray(sample_times, dtype=float))
    out = np.empty((st.size, 2), dtype=float)
    t_max = float(st[-1])
    M, P = float(M0), float(P0)
    t = 0.0
    si = 0
    while t < t_max:
        fM = 1.0 if kappa_M == 0.0 else 1.0 / (1.0 + (kappa_M * P) ** n_M)
        gP = 1.0 if kappa_P == 0.0 else 1.0 / (1.0 + (kappa_P * P) ** n_P)
        a1 = alpha_M * (1.0 + S_a) * fM
        a2 = beta_M * M
        a3 = alpha_P * M * gP
        a4 = beta_P * P
        a0 = a1 + a2 + a3 + a4
        if a0 <= 0.0:
            break
        t_next = t + rng.exponential(1.0 / a0)
        while si < st.size and st[si] <= t_next:
            out[si] = (M, P)
            si += 1
        t = t_next
        r = rng.random() * a0
        if r < a1:
            M += 1.0
        elif r < a1 + a2:
            M -= 1.0
        elif r < a1 + a2 + a3:
            P += 1.0
        else:
            P -= 1.0
    while si < st.size:
        out[si] = (M, P)
        si += 1
    return out


def ensemble_cv(*, n_runs, sample_times, seed, M0, P0, **params):
    """Run an ensemble and return means and coefficients of variation for mRNA and protein, pooled
    over all runs and sample times (after the burn-in implied by sample_times starting past 0)."""
    rng = np.random.default_rng(seed)
    Ms, Ps = [], []
    for _ in range(int(n_runs)):
        s = gillespie_nar(M0=M0, P0=P0, sample_times=sample_times, rng=rng, **params)
        Ms.append(s[:, 0])
        Ps.append(s[:, 1])
    M = np.concatenate(Ms)
    P = np.concatenate(Ps)
    return {
        "mean_M": float(M.mean()), "cv_M": float(M.std() / M.mean()),
        "mean_P": float(P.mean()), "cv_P": float(P.std() / P.mean()),
        "n_samples": int(M.size),
    }
