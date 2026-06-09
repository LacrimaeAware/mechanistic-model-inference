# Experiment 02: identifiability map for the three models

## Goal

For each model (M1 simplistic, M2 transcriptional, M3 translational) and each observation scheme
(protein only, or mRNA plus protein), determine which dimensional parameters are identifiable and name
the non-identifiable directions. Parameters: transcription rate k_m, mRNA decay d_m, translation rate
k_p, protein decay d_p, and the regulation strength kappa (kappa_M for M2, kappa_P for M3). The Hill
coefficient n = 2 and the input level S are held fixed (treated as known).

## Method

- Dimensional simulator: `make_dimensional_simulator` in `src/mechanistic_inference/models.py`, used
  through the carried-over identifiability tooling.
- Structural identifiability: Fisher information rank, and the eigenvector of the smallest eigenvalue
  (the non-identifiable direction). The rank verdict does not depend on the noise level.
- Practical identifiability: simulate normalized noisy data (3% relative noise, 80 time points), fit by
  maximum likelihood, profile each parameter, and classify by whether the profile rises above the 95%
  threshold for one parameter (delta = 1.92) on both sides.
- Script: `experiments/02_identifiability/identifiability_map.py` (writes
  `results/fig04_profile_likelihood.png` and `results/identifiability_map_summary.txt`). Tests:
  `tests/test_identifiability_nar.py`.

## Results

Structural map (Fisher rank, true parameters from Table 1, 1/kappa = 50 molecules):

| model | protein only | mRNA + protein | non-identifiable direction (protein only) |
| ----- | ------------ | -------------- | ----------------------------------------- |
| M1    | rank 3 / 4   | rank 4 / 4     | k_m vs k_p (transcription/translation product) |
| M2    | rank 4 / 5   | rank 5 / 5     | k_m vs k_p (kappa_M identifiable)         |
| M3    | rank 4 / 5   | rank 5 / 5     | k_m vs k_p (kappa_P identifiable)         |

The single non-identifiable direction under protein-only observation is the same in all three models:
the transcription/translation rate product. Observing mRNA restores full rank in every case. The
regulation strength kappa is identifiable from protein alone.

Analytic reason for the degeneracy: scaling k_m -> k_m c and k_p -> k_p / c rescales the mRNA
trajectory by c and leaves the protein equation (k_p M) unchanged, so the protein trajectory is
identical. This holds for the regulated models too, because the scaling does not touch the Hill term.
The test `test_rate_product_degeneracy_leaves_protein_invariant` checks this directly.

Practical map (M3, 3% noise): k_m is flat (non-identifiable) under protein-only and bounded once mRNA
is observed; kappa_P is bounded under both schemes, but the profile is shallow, so its confidence
interval is wide. The structural verdict (kappa identifiable) and the practical verdict (kappa
identifiable but weakly constrained) agree in direction and differ in strength.

## Hypotheses, with evidence for and against

- Protein-only leaves k_m and k_p non-identifiable in M1. Evidence for: matches the carried-over
  baseline and the analytic scaling symmetry. Evidence against: none found. Status: confirmed.
- Transcriptional feedback (M2) might break the product degeneracy under protein-only by feeding the
  protein back into mRNA production. Evidence for: f(P) makes mRNA production depend on P and kappa_M.
  Evidence against (decisive): the scaling k_m c, k_p / c leaves the protein trajectory invariant
  regardless of f(P), and the Fisher null direction is still the k_m vs k_p product. Status: refuted.
- kappa is identifiable from protein alone. Evidence for: full Fisher rank in that coordinate and a
  bounded profile. Caveat: the kappa_P profile is shallow at 3% noise and 80 samples, so the interval
  is wide. Status: supported structurally, weak practically.

## Methodology audit

The verdicts were checked five ways against the risk of method artifacts
(`experiments/02_identifiability/verify_methodology.py`, output in
`results/verify_methodology_summary.txt`).

1. Spectral gap: the non-identifiable eigenvalue is 1e-16 or exactly 0 relative to the largest, while
   the smallest identifiable eigenvalue is 1e-5 to 1e-4 relative. The rank is a clear call, not a
   threshold-sensitive one.
2. Finite-difference robustness: the M1 protein-only null direction aligns with the analytic k_m/k_p
   direction ([+1, 0, -1, 0]/sqrt(2)) to |dot| = 1.0000 for step sizes 1e-4, 1e-5, 1e-6. The null
   direction is the analytic degeneracy, not numerical noise.
3. Method agreement: for M3 with both channels (identifiable), the Fisher standard error and the
   profile-likelihood interval agree in order and direction. The profile half-width at the 95%
   threshold is about 1.3 to 1.5 times the Gaussian prediction (1.96 x Fisher SE), i.e. slightly wider,
   as expected when the likelihood has non-Gaussian tails. The profile half-widths are quantized by the
   grid resolution (0.2 in log units), so the comparison is order-of-magnitude, not exact.
4. MLE recovery: fitting both channels recovers all five parameters (max log error 0.07). Fitting
   protein only, started off the truth along the degenerate direction, recovers the k_m k_p product
   (log error 0.06) but not the split (log k_m error 0.50), and still recovers kappa_P (error 0.08).
   The degeneracy is reproduced operationally, not just asserted.
5. Noise sweep (M3): k_m is non-identifiable at every noise level tested (0.5% to 30%), confirming the
   degeneracy is structural, not noise-driven. kappa_P crosses from identifiable to non-identifiable
   between 1% and 3% noise under protein-only, and between 3% and 10% with mRNA observed, so the mRNA
   channel buys roughly a decade of noise tolerance for kappa. The crossover is approximate (one noise
   realization per setting, coarse sigma grid).

Net: the structural conclusion is corroborated by three independent routes (the analytic scaling
symmetry, the Fisher rank, and MLE recovery). The practical conclusion is consistent across the Fisher
standard error, the profile width, and the noise sweep, with the quantitative thresholds approximate.

## Open items

- Sampling design: the noise sweep is done (above), but sampling density and horizon were held fixed.
  Open hypothesis, untested: an input step that drives the protein through the Hill midpoint tightens
  the kappa interval.
- M2 practical profiles (only M3 was profiled in detail).
- The single-sigma likelihood weights mRNA and protein equally after per-channel normalization; a more
  realistic per-channel noise model is a later refinement.

## Notes

- The structural conclusion (which directions are null) is independent of the noise level; the
  practical conclusion is not.
- This experiment answers parameter identifiability within each model. Telling the models apart from
  each other (M1 vs M3, M2 vs M3) is the separate model-discrimination question (Phase 4).
