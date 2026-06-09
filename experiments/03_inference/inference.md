# Experiment 03: inference under noise (Phase 3)

## Goal

Move from "which parameters are identifiable" to "how well are they recovered from noisy data, and what
does the uncertainty look like." Two parts: repeated maximum-likelihood recovery as a function of
noise, and the shape of the Bayesian posterior.

## Method

- Recovery (M3, dimensional parameters): for each noise level and observation scheme, fit 12 noisy
  realizations from a perturbed start and report the median absolute log-error per parameter, plus the
  k_m k_p product error and the k_m split error.
- Posterior (M1, MCMC): sample the posterior with emcee (32 walkers, 4000 steps, 1000 discarded) using
  the closed-form cascade so the sampling is fast and exact. Flat prior in log space, Gaussian
  likelihood with known sigma. Report the correlation of log k_m and log k_p and the spread along the
  product direction (sum of logs) versus the ridge direction (difference of logs).
- Script: `experiments/03_inference/recovery_and_posterior.py` (writes `results/fig05_posterior_ridge.png`
  and `results/inference_summary.txt`).

## Results

Recovery (M3), median absolute log-error:

| scheme        | noise | k_m   | k_p   | kappa_P | k_m*k_p product | k_m split |
| ------------- | ----- | ----- | ----- | ------- | --------------- | --------- |
| protein-only  | 1%    | 0.216 | 0.141 | 0.084   | 0.030           | 0.216     |
| protein-only  | 5%    | 0.202 | 0.413 | 0.238   | 0.245           | 0.202     |
| mRNA+protein  | 1%    | 0.017 | 0.044 | 0.042   | 0.038           | 0.017     |
| mRNA+protein  | 5%    | 0.089 | 0.187 | 0.216   | 0.216           | 0.089     |

With both channels, every parameter is recovered (k_m to 0.017 log units at 1% noise). With protein
only, the k_m k_p product is recovered (0.030) while k_m itself is not pinned (0.216, about ten times
worse). Caveat: under protein-only the split error largely reflects the optimizer staying near its
perturbed starting point along the flat ridge, so this number is start-dependent and understates the
true freedom; the posterior below is the cleaner demonstration.

Posterior geometry (M1):

| scheme        | corr(log k_m, log k_p) | std(sum) | std(diff) | ratio diff/sum |
| ------------- | ---------------------- | -------- | --------- | -------------- |
| protein-only  | -0.93                  | 0.80     | 4.24      | 5.3            |
| mRNA+protein  | -0.26                  | 0.074    | 0.084     | 1.1            |

Under protein only the posterior is a long ridge along log k_m + log k_p = const: the product is
constrained (spread 0.80) and the split is not (spread 4.24, hitting the prior bounds). Observing mRNA
closes it into a tight, nearly round blob (figure `fig05_posterior_ridge.png`). This is the same
degeneracy seen by Fisher rank and MLE recovery, now in the full posterior.

## Hypotheses, with evidence for and against

- The protein-only degeneracy appears as a posterior ridge. Evidence for: correlation -0.93, ridge
  ratio 5.3, and the product spread (0.80) much smaller than the split spread (4.24). Evidence against:
  none. Status: confirmed, and it agrees with the Fisher and MLE results (a fourth independent route).
- Observing mRNA gives bounded recovery of all parameters. Evidence for: full-rank posterior blob
  (ratio 1.1, low correlation) and recovery to ~0.02-0.04 log units at 1% noise. Status: confirmed.

## Open items

- The recovery split error under protein-only is start-dependent (optimizer on a flat ridge); a cleaner
  metric is the posterior spread, already reported. Worth replacing the point-estimate metric with the
  posterior in a later pass.
- MCMC was run on M1 (closed form) for speed. The regulated models should show the same ridge plus a
  bounded kappa direction; running MCMC on M3 (ODE) is the natural confirmation, deferred for runtime.
- Convergence was not formally diagnosed (no R-hat / autocorrelation report); the ridge-vs-blob
  conclusion is qualitative and robust, but a quantitative posterior would need convergence checks.

## Notes

- "Identifiable" continues to mean a finite confidence region, not a tight one; the protein-only product
  is constrained but still has a finite width that grows with noise.
