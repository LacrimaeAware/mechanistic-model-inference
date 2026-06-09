# Methodology review

A higher-level check than data verification: for each phase, is the method the right tool, applied
correctly, with its assumptions met? This is about the approach, not re-running numbers. Findings on
specific results live in `docs/verification_debt.md`; this is about whether the methods themselves are
sound. External corroboration of "standard methods" is pending the literature research.

## Phase 1 - reproduction (forward simulation)

Method: integrate the ODEs, compare to the paper's reported behavior and to an analytic solution.
Appropriate and unambiguous. Validated by an analytic cross-check (M1 matches its closed form to 3e-8),
which confirms the solver. Assumption met: the equations were transcribed from the source and confirmed
against the rendered PDF. Limitation: the comparison to the paper's figures is qualitative (shapes and
orderings), not digitized point-by-point. Verdict: sound.

## Phase 2 - identifiability (Fisher rank + profile likelihood) - the strongest part

Method: structural identifiability via Fisher information rank; practical identifiability via profile
likelihood. These are the standard tools for this question. Applied correctly: the Fisher matrix is a
finite-difference Jacobian (step-robustness checked), the rank uses a relative-eigenvalue threshold with
a verified clear spectral gap, and the profile re-optimizes the other parameters. Crucially, the
structural result is also proven analytically (the k_m c, k_p/c scaling leaves the protein trajectory
invariant), so the conclusion does not rest on the numerical method at all. Assumptions: Gaussian noise;
Fisher is local, but the analytic proof makes the result global. Verdict: sound and method-independent -
the most trustworthy part of the project.

## Phase 3 - inference (MLE recovery + MCMC)

Method: maximum-likelihood recovery over noise levels; MCMC for the posterior. Standard. The MLE recovery
agrees with Phase 2 (product recovered, split not, from protein). Pitfall encountered and caught: the M3
MCMC does not converge on the flat ridge (steps/tau below threshold), a known hard case for ensemble
samplers. Verdict: sound for M1; the M3 posterior is an acknowledged methodological gap (needs nested
sampling or analytic marginalization), not a result.

## Phase 4 - discrimination (AIC selection rate) - a methodology error, since corrected

Method as first applied: fit two candidate models, compare AIC, repeat over noise draws, report the
selection rate. This was the WRONG tool for the question. When two non-nested models fit the data nearly
identically (which a deterministic check showed: mutual-fit residual under 1%, below the 3% noise), the
AIC winner is decided by sub-noise wobble, so the selection rate is an unstable quantity with no stable
value. The symptom was the rate swinging between runs; the cause was applying a discrimination statistic
where the models are structurally near-indistinguishable. The correct method is the deterministic
structural check (can one model reproduce another's noise-free output, and by how much relative to the
measurement noise). Verdict: the original AIC-rate method is retracted for the protein-only cases; the
structural statement (protein cannot separate the mechanisms, mRNA can) stands. This is the project's
main methodological lesson.

## Phase 5 - real-data fit (least squares + Fisher rank on real data)

Method: fit the closed-form translation model to real GFP data; check identifiability by Fisher rank.
Appropriate. Externally corroborated: the Fisher rank recovers the published scale*k*m0 non-identifiability
(Pieschner et al. 2022), which is independent validation that the method is applied correctly on real
data. Limitations: mean trace (not single-cell mixed effects), single observation channel, the delta/gamma
label symmetry, one noise model. Verdict: sound for the structural claim; the parameter estimates are a
first pass.

## Phase 7 - experimental design (deterministic Fisher)

Method: Fisher information at the true parameters to score observation designs. Deterministic, so stable
(the antidote to the Phase 4 instability). Builds directly on Phase 2. Limitation: finite-difference
Fisher (accuracy), and full rank is structural - the practical error with one mRNA point is not yet
quantified. Verdict: sound.

## Overall

The identifiability core (Phases 2, 7) and the reproduction (Phase 1) use appropriate, standard,
correctly-applied methods, with the strongest result backed by an analytic proof rather than a numerical
method. The real-data fit (Phase 5) is method-validated against a published result. The two methodological
weaknesses are both documented: the M3 MCMC convergence (Phase 3) and the AIC-discrimination instability
(Phase 4, corrected to a structural statement). The recurring failure mode to guard against is applying a
noisy statistic to a near-degenerate problem and reading its fluctuations as results.
