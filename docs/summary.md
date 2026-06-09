# Project summary: the story so far

Built from the result files in `results/` (the actual numbers, not the prose). One through-line connects
every phase.

## The one finding everything points to

For these mRNA -> protein autoregulation models, the protein channel alone is information-poor; the mRNA
channel carries the information. Every phase says this from a different angle, and the practical message
is the same: to pin down the parameters or the mechanism, measure mRNA.

## Phase by phase, from the actual numbers

**Phase 1 (reproduction).** The three models reproduce the source paper: only the transcriptional model
oscillates (oscillation discriminant -12.7 vs +0.16; M1/M3 overshoot ~1e-13), autoregulation speeds the
response (half-response time 13.2 for M1 vs 3.6-4.1 for M2/M3), and the stochastic protein-noise ordering
is M2 > M1 > M3 (CV 0.153 > 0.141 > 0.091; paper 0.16/0.14/0.10). Two structural facts fall out and drive
everything after: **M3 leaves mRNA identical to M1** (translational regulation does not touch mRNA; max
difference 6e-9), and **M2 and M3 share the same protein steady state** under symmetric feedback (both
1.333). [solid: analytic + tests]

**Phase 2 (identifiability) - the solid core.** From protein alone, the transcription rate and the
translation rate are not separately identifiable in any of the three models (Fisher rank 3/4 for M1, 4/5
for M2/M3; the null direction is the k_m vs k_p product). The regulation strength kappa is structurally
identifiable but only marginally so from protein at realistic noise (identifiable in 5/10 and 6/10 noise
draws). Observing mRNA restores full rank. Verified five independent ways: Fisher rank, an analytic
scaling-symmetry proof, MLE recovery (product error 0.06 but split error 0.50 from protein), a
step-robustness check, and a noise sweep (k_m non-identifiable at every noise level; kappa degrades from
identifiable to not as noise grows). [solid]

**Phase 3 (inference) - confirms Phase 2.** Fitting noisy data recovers the k_m*k_p product (error 0.03
at 1% noise) but not the split (0.22) from protein; both channels recover everything (0.02). The M1
posterior is a literal ridge under protein-only (correlation -0.93, long axis 5x the short) that mRNA
closes into a blob (correlation -0.26). The M3 posterior is the same shape but its sampler did not
converge (reported, not used). [solid for M1; M3 MCMC under-converged]

**Phase 4 (discrimination) - read correctly, the same story.** The discriminating information lives in
the mRNA channel because that is where the models differ structurally: M3's mRNA equals M1's, so mRNA
cannot separate M1 from M3 (M1 wins by the parsimony penalty, delta-AIC +2.0); M2's mRNA is strongly
suppressed (steady-state mRNA ~1.0 vs ~3.0 for M1/M3), so mRNA separates M2 from M3 decisively (delta-AIC
~1000). Protein alone does NOT reliably discriminate: a deterministic check shows any model re-fits
another's clean protein curve to under 1% RMS (M2->M3 0.055%, M1->M3 0.87%), far below the 3% measurement
noise, so the protein-only AIC winner is decided by noise. The Monte-Carlo protein-only selection rates I
reported (numbers from 16% to 88% across runs) were measuring this sub-noise tie and are retracted; the
stable statement is structural - protein cannot separate the mechanisms, mRNA can. [mRNA-driven results
solid; protein-only rates retracted]

**Phase 5 (real data, provisional).** On a real single-cell mRNA->protein dataset (Frohlich 2018, open-
loop translation), the translation model fits the mean trace to 1.4% and the pipeline recovers the
published structural non-identifiability (the scale*k*m0 product; Fisher rank 5/7 with both null
directions in that block). A first pass, not independently re-verified, and not the source's
autoregulation system. [provisional]

## Solid vs not

- Solid: Phase 1 reproduction; Phase 2 identifiability; Phase 3 for M1; the structural discrimination
  facts (mRNA separates the models, protein does not, with the deterministic residuals to back it).
- Provisional or retracted: all Monte-Carlo discrimination selection rates (sub-noise, unstable); the M3
  MCMC (under-converged); the real-data fit (one pass). See `docs/verification_debt.md`.

## Takeaway

One consistent, defensible message across every phase: for these systems, protein measurements alone
cannot pin down the rate parameters or the regulatory mechanism; the mRNA channel is what makes both
identifiable and the mechanism distinguishable. That is the experimental-design conclusion and the spine
a write-up would be built on.
