"""Mechanistic ODE models for the project's first target.

Source (transcribed from the paper text and confirmed against the rendered PDF, not from memory;
full citation in docs/literature.md): Ryzowicz and Yildirim (2023), "Differential roles of
transcriptional and translational negative autoregulations in protein dynamics." An input signal
up-regulates a protein's transcription; the protein then represses its own production either
transcriptionally or translationally.

Three nested models are defined. In dimensional form (Eqs. 1-3) they share

    dM/dt = alpha_M f(P) - beta_M M ,    dP/dt = alpha_P M g(P) - beta_P P ,
    f(P) = 1 / (1 + (kappa_M P)^n_M) ,   g(P) = 1 / (1 + (kappa_P P)^n_P) .

Nondimensionalizing with tau = beta_M t, m = (beta_M/alpha_M) M, p = (beta_M^2/(alpha_M alpha_P)) P,
theta = beta_P/beta_M, and kappa_tilde = (alpha_M alpha_P / beta_M^2) kappa (so kappa_tilde p = kappa P
exactly) gives the dimensionless system implemented here (Eqs. 4-7):

    dm/dtau = (1 + S_a) F(p) - m ,       dp/dtau = m G(p) - theta p ,
    F(p) = 1 / (1 + (kappa_tilde_M p)^n_M) ,   G(p) = 1 / (1 + (kappa_tilde_P p)^n_P) .

The three models are the parameter limits:
    M1 simplistic        : kappa_tilde_M = kappa_tilde_P = 0  (F = G = 1)
    M2 transcriptional   : kappa_tilde_M > 0, kappa_tilde_P = 0  (repression on transcription)
    M3 translational     : kappa_tilde_M = 0, kappa_tilde_P > 0  (repression on translation)

Two ambiguities in the linearized text layer were resolved by reading the rendered PDF pages:
- The Hill term is (kappa P) raised to n, i.e. 1/(1 + (kappa P)^n), not kappa P^n (Eq. 3, p. 4).
- The input signal scales the transcription rate multiplicatively: the baseline dimensionless
  transcription rate 1 becomes (1 + S_a), so the transcription term is (1 + S_a) F(p). S_a in [1, 4]
  corresponds to a 2- to 5-fold increase, matching the paper's text (p. 15-16). For M1 and M3 this
  reduces to (1 + S_a) since F = 1.

Table-1 dimensional parameters: alpha_M = 1.04, alpha_P = 7.70 /min, beta_M = 0.35 /min (2-min mRNA
half-life), beta_P = 0.02 /min (30-min protein half-life), n_M = n_P = 2, inhibition strengths
1/kappa in [10, 500] molecules, S_a in [1, 4]. These give theta = beta_P/beta_M ~= 0.0571; the
oscillation analysis (Figs. 3-4) instead sweeps theta in {0.2, 0.4, 0.6}.
"""
from __future__ import annotations

import numpy as np


def _rhs(tau, y, theta, kM, kP, nM, nP, Sa):
    """Right-hand side of the dimensionless system (Eqs. 4-7) for the general model."""
    m, p = y
    F = 1.0 if kM == 0.0 else 1.0 / (1.0 + (kM * p) ** nM)
    G = 1.0 if kP == 0.0 else 1.0 / (1.0 + (kP * p) ** nP)
    return ((1.0 + Sa) * F - m, m * G - theta * p)


def simulate_nar(t, *, theta, kappa_tilde_M=0.0, kappa_tilde_P=0.0, n_M=2.0, n_P=2.0, S_a=0.0,
                 m0=0.0, p0=0.0, rtol=1e-9, atol=1e-12):
    """Integrate the dimensionless model. Returns an array of shape (len(t), 2) = columns [m, p].

    The model is selected by the inhibition strengths: kappa_tilde_M/P = 0 turns off the
    corresponding Hill regulation. ``theta`` and the kappa_tilde values are already dimensionless.
    The system is mildly stiff (fast mRNA, slow protein when theta is small), so LSODA is used to
    match the paper's stiff solver (ode15s).
    """
    from scipy.integrate import solve_ivp

    t = np.asarray(t, dtype=float)
    sol = solve_ivp(_rhs, (float(t[0]), float(t[-1])), [m0, p0], t_eval=t,
                    args=(theta, kappa_tilde_M, kappa_tilde_P, n_M, n_P, S_a),
                    method="LSODA", rtol=rtol, atol=atol)
    return sol.y.T


def simulate_simplistic(t, *, theta, S_a=0.0, **kwargs):
    """Model M1 (i): no regulation. dm/dtau = (1 + S_a) - m ; dp/dtau = m - theta p."""
    return simulate_nar(t, theta=theta, S_a=S_a, **kwargs)


def simulate_transcriptional_nar(t, *, theta, kappa_tilde_M, n_M=2.0, S_a=0.0, **kwargs):
    """Model M2 (ii): transcriptional negative autoregulation (repression of mRNA production)."""
    return simulate_nar(t, theta=theta, kappa_tilde_M=kappa_tilde_M, n_M=n_M, S_a=S_a, **kwargs)


def simulate_translational_nar(t, *, theta, kappa_tilde_P, n_P=2.0, S_a=0.0, **kwargs):
    """Model M3 (iii): translational negative autoregulation (repression of protein production)."""
    return simulate_nar(t, theta=theta, kappa_tilde_P=kappa_tilde_P, n_P=n_P, S_a=S_a, **kwargs)


def simplistic_analytic(t, *, theta, S_a=0.0, m0=0.0, p0=0.0):
    """Closed-form solution of M1 (Eq. 8), generalized to a nonzero signal (mRNA steady state 1+S_a).

    Requires theta != 1 (the theta == 1 case has a different, degenerate form). Returns shape
    (len(t), 2). Used as the analytic ground truth for the numerical solver.
    """
    if theta == 1.0:
        raise ValueError("Eq. 8 closed form requires theta != 1.")
    t = np.asarray(t, dtype=float)
    m_ss = 1.0 + S_a
    m = m_ss + (m0 - m_ss) * np.exp(-t)
    A = (m0 - m_ss) / (theta - 1.0)
    p = A * np.exp(-t) + (p0 - A - m_ss / theta) * np.exp(-theta * t) + m_ss / theta
    return np.column_stack([m, p])


def steady_state(*, theta, kappa_tilde_M=0.0, kappa_tilde_P=0.0, n_M=2.0, n_P=2.0, S_a=0.0):
    """Steady state (m*, p*) of the dimensionless model.

    M1 and M3 have m* = 1 + S_a analytically; M2 reduces transcription, lowering m*. In every case
    p* is the unique root of the scalar steady-state equation, found from the M1 guess p = (1+S_a)/theta.
    """
    from scipy.optimize import brentq

    def F(p):
        return 1.0 if kappa_tilde_M == 0.0 else 1.0 / (1.0 + (kappa_tilde_M * p) ** n_M)

    def G(p):
        return 1.0 if kappa_tilde_P == 0.0 else 1.0 / (1.0 + (kappa_tilde_P * p) ** n_P)

    # At steady state m* = (1 + S_a) F(p*) and m* G(p*) = theta p*, so
    # (1 + S_a) F(p*) G(p*) - theta p* = 0 in the single unknown p*.
    def resid(p):
        return (1.0 + S_a) * F(p) * G(p) - theta * p

    hi = (1.0 + S_a) / theta  # p* for the unregulated case; regulation only lowers p*, so this brackets
    p_star = brentq(resid, 1e-12, hi + 1e-9)
    m_star = (1.0 + S_a) * F(p_star)
    return m_star, p_star


def _translation_signal(amplitude, delta, gamma, t0, offset, t):
    """Closed-form GFP signal of the mRNA-bolus translation model (Frohlich/Leonhardt), for experiment 05.

    For t >= t0 (tau = t - t0): GFP(tau) = amplitude/(gamma - delta) * (exp(-delta tau) - exp(-gamma tau)),
    observed as amplitude*... + offset; before t0 the signal is the baseline offset. ``amplitude`` is the
    lumped scale*k*m0 (see the two simulators below).
    """
    t = np.asarray(t, dtype=float)
    tau = np.maximum(t - t0, 0.0)
    if abs(gamma - delta) < 1e-9:
        g = tau * np.exp(-delta * tau)                     # limit as gamma -> delta
    else:
        g = (np.exp(-delta * tau) - np.exp(-gamma * tau)) / (gamma - delta)
    y = amplitude * g + offset
    return np.where(t < t0, offset, y)


def simulate_translation(theta, t):
    """Translation model in the IDENTIFIABLE parameterization. theta = log[amplitude, delta, gamma, t0,
    offset]; returns shape (len(t), 1) (GFP is the single observable). ``amplitude`` is the lumped
    scale*k*m0 (only this product is identifiable from GFP alone)."""
    amplitude, delta, gamma, t0, offset = np.exp(np.asarray(theta, dtype=float))
    return _translation_signal(amplitude, delta, gamma, t0, offset, t).reshape(-1, 1)


def simulate_translation_full(theta, t):
    """Translation model in the SEPARATED parameterization, for the identifiability demonstration.
    theta = log[scale, k, m0, delta, gamma, t0, offset]; amplitude = scale*k*m0. The three factors scale,
    k, m0 enter only through their product, so two of their directions are non-identifiable from GFP
    alone (the published result). Returns shape (len(t), 1)."""
    scale, k, m0, delta, gamma, t0, offset = np.exp(np.asarray(theta, dtype=float))
    return _translation_signal(scale * k * m0, delta, gamma, t0, offset, t).reshape(-1, 1)


def make_dimensional_simulator(model="M1", *, n=2.0, S=0.0, m0=0.0, p0=0.0,
                               rtol=1e-8, atol=1e-10):
    """Build a ``simulate(theta, t) -> (len(t), 2)`` closure for the identifiability tooling.

    The closure integrates the dimensional model in molecule units, so the transcription/translation
    rates are explicit and the protein-only product degeneracy (the baseline result) can appear:

        dM/dt = k_m (1 + S) f(P) - d_m M ,   dP/dt = k_p M g(P) - d_p P .

    ``theta`` is the log of the free parameters, in a fixed order per model:
        M1 (no regulation)        : [k_m, d_m, k_p, d_p]
        M2 (transcriptional NAR)  : [k_m, d_m, k_p, d_p, kappa_M]   (f(P) = 1/(1+(kappa_M P)^n))
        M3 (translational NAR)    : [k_m, d_m, k_p, d_p, kappa_P]   (g(P) = 1/(1+(kappa_P P)^n))
    The Hill coefficient ``n`` and the input level ``S`` are held fixed (treated as known), matching the
    paper's fixed n = 2; identifiability is asked about the rate and regulation parameters.
    """
    from scipy.integrate import solve_ivp

    model = model.upper()
    if model not in ("M1", "M2", "M3"):
        raise ValueError("model must be 'M1', 'M2', or 'M3'")

    def simulate(theta, t):
        p = np.exp(np.asarray(theta, dtype=float))
        k_m, d_m, k_p, d_p = p[0], p[1], p[2], p[3]
        kappa_M = p[4] if model == "M2" else 0.0
        kappa_P = p[4] if model == "M3" else 0.0

        def rhs(_t, y):
            M, P = y
            f = 1.0 if kappa_M == 0.0 else 1.0 / (1.0 + (kappa_M * P) ** n)
            g = 1.0 if kappa_P == 0.0 else 1.0 / (1.0 + (kappa_P * P) ** n)
            return (k_m * (1.0 + S) * f - d_m * M, k_p * M * g - d_p * P)

        tt = np.asarray(t, dtype=float)
        try:
            sol = solve_ivp(rhs, (float(tt[0]), float(tt[-1])), [m0, p0], t_eval=tt,
                            method="LSODA", rtol=rtol, atol=atol)
            y = np.asarray(sol.y, dtype=float)
            if sol.success and y.shape == (2, tt.size):
                return y.T
        except Exception:
            pass
        # Integration failed (extreme parameters during optimization): return a large finite sentinel
        # so the likelihood is huge-but-finite and the optimizer steers away from this region.
        return np.full((tt.size, 2), 1e6, dtype=float)

    return simulate


def m2_oscillation_discriminant(*, theta, kappa_tilde_M, n_M=2.0, S_a=0.0):
    """Discriminant of the M2 Jacobian eigenvalues (Eq. 9): Delta = (theta - 1)^2 - 4 h1(p*).

    Delta < 0 means complex eigenvalues, i.e. damped oscillation; Delta >= 0 means a monotone
    (overdamped) approach. h1(p*) = n_M (kappa_tilde_M p*)^n_M / [ (1 + (kappa_tilde_M p*)^n_M)^2 p* ].
    """
    _, p_star = steady_state(theta=theta, kappa_tilde_M=kappa_tilde_M, n_M=n_M, S_a=S_a)
    x = (kappa_tilde_M * p_star) ** n_M
    h1 = n_M * x / ((1.0 + x) ** 2 * p_star)
    return (theta - 1.0) ** 2 - 4.0 * h1
