# Experiment 09: structural distinguishability of transcriptional vs translational autoregulation

## Goal

Answer the one genuinely open question the literature research identified (the source paper builds M2/M3
but never asks it): can transcriptional (M2) and translational (M3) negative autoregulation be
distinguished from time-course data, and from which channel? Answered with a deterministic method so the
result is stable - the correction of the earlier noisy AIC selection rates, which were retracted.

## Method

For each feedback strength and observation channel, fit each model to the OTHER model's NOISE-FREE output
by best-effort multi-start least squares, and record the leftover RMS residual as a percent of the peak.
That residual is the structural distance between the mechanisms after the wrong model does its best to
mimic the right one. If the residual is below the measurement noise, the mechanisms are indistinguishable
from that channel; if above, they are distinguishable. No noise draws are involved, so the numbers do not
change between runs. Script: `experiments/09_structural_distinguishability/structural_distinguishability.py`,
figure `results/fig13_structural_distinguishability.png`.

## Results (deterministic)

Best wrong-mechanism mimic residual (% of peak); "min" is the harder-to-distinguish direction:

| channel       | 1/kappa | M2 mimics M3 | M3 mimics M2 | min   |
| ------------- | ------- | ------------ | ------------ | ----- |
| protein-only  | 200     | 0.060%       | 0.126%       | 0.060%|
| protein-only  | 50      | 0.055%       | 0.290%       | 0.055%|
| protein-only  | 20      | 0.025%       | 0.236%       | 0.025%|
| protein-only  | 10      | 0.037%       | 0.250%       | 0.037%|
| protein-only  | 5       | 0.040%       | 1.358%       | 0.040%|
| mRNA+protein  | 200     | 0.624%       | 8.455%       | 0.624%|
| mRNA+protein  | 50      | 1.348%       | 8.190%       | 1.348%|
| mRNA+protein  | 20      | 1.718%       | 7.279%       | 1.718%|
| mRNA+protein  | 10      | 1.933%       | 6.808%       | 1.933%|
| mRNA+protein  | 5       | 2.084%       | 6.685%       | 2.084%|

- Protein-only: the best mimic residual is 0.03-0.06% of peak at every feedback strength, far below any
  realistic noise (1-10%). The mechanisms are structurally indistinguishable from protein.
- mRNA+protein: asymmetric. M3 cannot mimic M2 (residual 6.7-8.5%, well above noise), so transcriptional
  autoregulation is clearly detectable with mRNA - only M2 regulates mRNA. But M2 CAN mimic M3 to 0.6-2%
  (the "min" column), so translational autoregulation is only marginally distinguishable even with mRNA,
  near the noise floor; it becomes easier as feedback strengthens (0.6% at 1/kappa=200 to 2.1% at 5).

## Reading (experimental design)

To tell which autoregulation mechanism a gene uses: protein measurement cannot do it at all. Measuring
mRNA reveals transcriptional autoregulation cleanly, but confirming translational autoregulation requires
low measurement noise (below ~1-2%) and is helped by strong feedback. This is the actionable answer to
the open question, stated structurally (a residual relative to noise), not as an unstable selection rate.

## Verification targets

- Whether the multi-start reaches the true global-best mimic at each cell (a worse local optimum would
  OVER-state distinguishability); more starts and a second optimizer would confirm the residuals.
- The percent-of-peak residual uses per-channel normalization; an absolute per-channel noise model would
  rescale the comparison to the noise lines.
- A rigorous structural-identifiability tool (SIAN / StructuralIdentifiability.jl) on M2 and M3 would
  independently confirm the protein-only structural indistinguishability.
