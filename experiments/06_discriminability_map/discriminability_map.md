# Experiment 06: discriminability map (transcriptional vs translational autoregulation)

## Goal

Answer the question the source paper raises by comparing the two mechanisms but does not address: across
feedback strength and observation channel, how often can data correctly identify which mechanism
(transcriptional M2 or translational M3) generated it? This is the source-relevant, experimental-design
result the project is aiming for.

## Method

For each cell (data-generating model, feedback strength, observation channel), generate noisy data over
N draws, fit both candidate models by multi-start least squares (over kappa, because the fit is stuck
near kappa = 0 otherwise), and record the selection rate: the fraction of draws the true model wins on
AIC. ~50% = indistinguishable; ~100% = reliably distinguishable. Script:
`experiments/06_discriminability_map/discriminability_map.py`; figure `results/fig10_discriminability_map.png`.

## A first pass that did not survive its own check (kept as a lesson)

The first run (N = 12) reported a striking result: under protein-only observation with the truth = M3
(translational), the true model won only 17-25% of the time, i.e. protein-only data appeared to point to
the WRONG mechanism. A direct check refuted this before it was reported:

- The M3 fit DOES reach the noise floor on M3 data (it is not failing), and M2 and M3 fit M3-generated
  protein data nearly identically (residual sum of squares tied to ~4 decimals). So the true protein-only
  rate should be near 50% (a coin flip between two near-equal fits), not 17%.
- The 17-25% was an artifact of two bugs: the RNG was seeded with `hash()` of a tuple, which Python
  randomizes per process (so the numbers were not even reproducible), and N = 12 gives huge variance on a
  ~50% rate. Both are fixed: reproducible integer seeds and N = 25.

The lesson, recorded because it will recur: a surprising number (here, a true model losing more than half
the time) is a prompt to check, not to report. See docs/verification_debt.md.

## Results (provisional; from the reproducible N = 25 re-run)

Robust, low-variance findings (rates at 0% or 100% have small sampling error and are mechanistically
sensible):

- mRNA + protein observation distinguishes the two mechanisms reliably (~100%) at every feedback
  strength. Mechanistically clear: only transcriptional autoregulation regulates mRNA, so the mRNA
  channel reveals the mechanism directly.
- The asymmetry is the real finding. Transcriptional autoregulation (truth M2) is mostly detectable from
  protein alone, more so as feedback strengthens. Translational autoregulation (truth M3) is close to
  indistinguishable from transcriptional using protein alone, because M2 and M3 fit M3 protein data
  almost identically.

Exact per-cell rates are in `results/discriminability_map_summary.txt`; the protein-only cells near 50%
are pinned only to about +/- 10% at N = 25 and should not be over-read.

Experimental-design reading (provisional): to identify which autoregulation mechanism a gene uses, the
mRNA channel is necessary; protein-only measurement can identify transcriptional autoregulation but not
translational, which it confounds with transcriptional.

## Verification targets (not yet done)

- Confirm the multi-start reaches the global optimum at each cell (more starts, a different optimizer).
- Cross-check the AIC selection against a likelihood-ratio or Bayes-factor calculation on a few cells.
- Sweep noise (only sigma = 0.03 here) and confirm the asymmetry is not an artifact of the fixed grid.
- Re-derive the "M2 mimics M3 in protein" claim analytically or with a second code path.
