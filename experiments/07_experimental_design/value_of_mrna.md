# Experiment 07: the value of the mRNA channel (deterministic experimental design)

## Goal

Turn the project's spine (protein alone is information-poor; mRNA carries the information) into a
quantitative, actionable experimental-design statement, using a deterministic metric (Fisher
information at the true parameters) rather than Monte-Carlo rates. Deterministic means stable: unlike the
discrimination selection rates, this result does not change between runs.

## Method

For each model, the Fisher information matrix is evaluated at the Table-1 parameters under each
observation scheme (per-channel normalized, sigma = 0.03 relative). Reported: rank, the smallest/largest
eigenvalue ratio, and the parameter standard errors (sqrt of the diagonal of the inverse Fisher matrix)
where the matrix is full rank. Then, with protein measured at all 80 timepoints, mRNA is added at K
timepoints (placed in the informative window, away from t = 0 where mRNA = 0) and the smallest K giving
full rank is found. Script: `experiments/07_experimental_design/value_of_mrna.py`.

## Results (deterministic)

- Protein-only is rank-deficient in every model (M1 3/4, M2 and M3 4/5). The deficient direction is the
  k_m vs k_p product (consistent with the Phase 2 analytic proof). One parameter combination therefore
  has infinite variance, and this is structural: no number of protein timepoints removes it.
- mRNA + protein is full rank in every model, with finite, modest standard errors (log units): M1
  k_m/d_m ~0.07, k_p/d_p ~0.02; M2 all ~0.03-0.07; M3 k_m/d_m ~0.07, k_p 0.16, d_p 0.11, kappa_P 0.14.
- A single informative mRNA timepoint restores full rank in all three models (the protein-only
  deficiency is rank one, so one independent mRNA measurement closes it). With evenly-spaced placement
  that includes t = 0 the experiment reports 2, because the t = 0 mRNA value is identically zero and
  carries no information; placed in the informative window, K = 1 suffices.

## Design conclusion (stable)

To estimate the rate parameters you cannot substitute dense protein sampling for mRNA: the protein-only
degeneracy is structural and survives any amount of protein data. A single well-timed mRNA measurement
removes it. The actionable statement is: measure mRNA at least once in the informative window; do not
spend the budget on more protein timepoints expecting to break the transcription/translation degeneracy.

## Verification targets (not yet done)

- The Fisher matrix is computed by finite differences; check against a tighter step and, where possible,
  the analytic sensitivities.
- "One timepoint suffices" should be checked across which timepoint is chosen (some placements are far
  more informative) and across parameter values, and the resulting standard error with a single mRNA
  point quantified (full rank is structural identifiability; the practical error with one point may still
  be large).
- The standard errors assume the per-channel normalized noise model; a per-channel absolute noise model
  would rescale them.
