"""Reproduce the deterministic behavior of the three negative-autoregulation models.

Source: Ryzowicz and Yildirim (2023); see docs/literature.md and src/mechanistic_inference/models.py.
Regenerates two figures into results/ and prints the numbers recorded in
experiments/01_reproduce_autoregulation_models/reproduce_autoregulation_models.md.

Run: python experiments/01_reproduce_autoregulation_models/reproduce.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from mechanistic_inference import models as M  # noqa: E402

RESULTS = os.path.join(os.path.dirname(__file__), "..", "..", "results")
os.makedirs(RESULTS, exist_ok=True)


def half_response_time(t, p):
    """Time for p to first reach halfway between p[0] and its final value."""
    target = p[0] + 0.5 * (p[-1] - p[0])
    above = np.where((p - target) * np.sign(p[-1] - p[0]) >= 0)[0]
    return float(t[above[0]]) if above.size else np.nan


def settling_time(t, p, frac=0.05):
    """Last time p leaves a +/- frac band around its final value (return/settling time)."""
    final = p[-1]
    band = frac * abs(final) if final != 0 else frac
    outside = np.where(np.abs(p - final) > band)[0]
    return float(t[outside[-1]]) if outside.size else 0.0


def main():
    lines = []

    def record(s):
        print(s)
        lines.append(s)

    # --- Correctness: M1 numerical vs analytic (Eq. 8) ---
    t = np.linspace(0, 50, 1000)
    num = M.simulate_simplistic(t, theta=0.0571, S_a=0.0, m0=0.0, p0=0.0)
    ana = M.simplistic_analytic(t, theta=0.0571, S_a=0.0, m0=0.0, p0=0.0)
    record(f"M1 numeric-vs-analytic (Eq. 8) max abs error: {np.max(np.abs(num - ana)):.2e}")

    # --- Figure 1: oscillation only in the transcriptional model (reproduces Fig. 4, theta=0.6) ---
    theta = 0.6
    t1 = np.linspace(0, 60, 3000)
    osc = M.simulate_transcriptional_nar(t1, theta=theta, kappa_tilde_M=2.0, n_M=8, S_a=0.0)  # 1/k=0.5
    mon = M.simulate_transcriptional_nar(t1, theta=theta, kappa_tilde_M=1 / 6, n_M=8, S_a=0.0)  # 1/k=6
    m3p = M.simulate_translational_nar(t1, theta=theta, kappa_tilde_P=2.0, n_P=8, S_a=0.0)
    d_osc = M.m2_oscillation_discriminant(theta=theta, kappa_tilde_M=2.0, n_M=8)
    d_mon = M.m2_oscillation_discriminant(theta=theta, kappa_tilde_M=1 / 6, n_M=8)
    record(f"M2 discriminant 1/kM=0.5: {d_osc:+.3f} (Delta<0 => damped oscillation)")
    record(f"M2 discriminant 1/kM=6.0: {d_mon:+.3f} (Delta>=0 => monotone)")
    record(f"M2 overshoot 1/kM=0.5: {osc[:, 1].max() - osc[-1, 1]:.4f}; "
           f"1/kM=3: {mon[:, 1].max() - mon[-1, 1]:.4f}; M3 (translational): {m3p[:, 1].max() - m3p[-1, 1]:.2e}")

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(t1, osc[:, 1], label=r"M2 transcriptional, $1/\tilde\kappa_M=0.5$ ($\Delta<0$)")
    ax.plot(t1, mon[:, 1], label=r"M2 transcriptional, $1/\tilde\kappa_M=6$ ($\Delta\geq0$)")
    ax.plot(t1, m3p[:, 1], "--", label=r"M3 translational, $1/\tilde\kappa_P=0.5$")
    ax.set_xlabel(r"dimensionless time $\tau$")
    ax.set_ylabel(r"dimensionless protein $p$")
    ax.set_title(r"Damped oscillation only in transcriptional autoregulation ($\theta=0.6$, $n=8$)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    f1 = os.path.join(RESULTS, "fig01_oscillation.png")
    fig.savefig(f1, dpi=140)
    record(f"wrote {os.path.relpath(f1)}")

    # --- Figure 2: faster response and quicker return under regulation (reproduces Fig. 5 ordering) ---
    theta = 0.0571
    kt = (1.04 * 7.70 / 0.35 ** 2) * (1.0 / 50.0)  # 1/kappa = 50 molecules -> dimensionless
    S_a = 3.0
    t2 = np.linspace(0, 400, 4000)
    # start each model at its own S_a = 0 steady state, then turn the signal on
    s1 = M.steady_state(theta=theta, S_a=0.0)
    s2 = M.steady_state(theta=theta, kappa_tilde_M=kt, n_M=2, S_a=0.0)
    s3 = M.steady_state(theta=theta, kappa_tilde_P=kt, n_P=2, S_a=0.0)
    p1 = M.simulate_simplistic(t2, theta=theta, S_a=S_a, m0=s1[0], p0=s1[1])[:, 1]
    p2 = M.simulate_transcriptional_nar(t2, theta=theta, kappa_tilde_M=kt, n_M=2, S_a=S_a, m0=s2[0], p0=s2[1])[:, 1]
    p3 = M.simulate_translational_nar(t2, theta=theta, kappa_tilde_P=kt, n_P=2, S_a=S_a, m0=s3[0], p0=s3[1])[:, 1]

    th = {"M1": half_response_time(t2, p1), "M2": half_response_time(t2, p2), "M3": half_response_time(t2, p3)}
    tr = {"M1": settling_time(t2, p1), "M2": settling_time(t2, p2), "M3": settling_time(t2, p3)}
    record(f"half-response time tau_h: M1={th['M1']:.1f}  M2={th['M2']:.1f}  M3={th['M3']:.1f}  "
           f"(M1 slowest: {th['M1'] > th['M2'] and th['M1'] > th['M3']})")
    record(f"return/settling time tau_r: M1={tr['M1']:.1f}  M2={tr['M2']:.1f}  M3={tr['M3']:.1f}  "
           f"(M1 slowest: {tr['M1'] > tr['M2'] and tr['M1'] > tr['M3']})")

    # normalize each protein trace to its own S_a=0 baseline (as the paper scales P by the baseline)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(t2, p1 / s1[1], label=f"M1 simplistic ($\\tau_h$={th['M1']:.0f})")
    ax.plot(t2, p2 / s2[1], label=f"M2 transcriptional ($\\tau_h$={th['M2']:.0f})")
    ax.plot(t2, p3 / s3[1], label=f"M3 translational ($\\tau_h$={th['M3']:.0f})")
    ax.set_xlabel(r"dimensionless time $\tau$")
    ax.set_ylabel(r"protein $p(\tau)\,/\,p(S_a{=}0)$")
    ax.set_title(r"Faster response under autoregulation (step $S_a=3$, $\theta=0.0571$)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    f2 = os.path.join(RESULTS, "fig02_response_speed.png")
    fig.savefig(f2, dpi=140)
    record(f"wrote {os.path.relpath(f2)}")

    # --- Identifiability-relevant degeneracies (for later phases) ---
    mrna_diff = np.max(np.abs(
        M.simulate_translational_nar(t2, theta=theta, kappa_tilde_P=kt, n_P=2, S_a=S_a, m0=s3[0], p0=s3[1])[:, 0]
        - M.simulate_simplistic(t2, theta=theta, S_a=S_a, m0=s1[0], p0=s1[1])[:, 0]))
    record(f"max|mRNA(M3) - mRNA(M1)| over the step: {mrna_diff:.2e} (mRNA cannot separate M1 from M3)")
    sym2 = M.steady_state(theta=0.3, kappa_tilde_M=1.5, n_M=2, S_a=1.0)[1]
    sym3 = M.steady_state(theta=0.3, kappa_tilde_P=1.5, n_P=2, S_a=1.0)[1]
    record(f"symmetric-parameter protein steady states: M2 p*={sym2:.6f}, M3 p*={sym3:.6f} "
           f"(coincide: {abs(sym2 - sym3) < 1e-6})")

    out = os.path.join(RESULTS, "reproduce_01_summary.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    record(f"wrote {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
