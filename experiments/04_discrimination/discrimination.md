# Experiment 04: model discrimination (Phase 4)

## Goal

Can the models be told apart from data, and under which measurement? This builds directly on the
fitting machinery: generate data from a known model, fit each candidate, and compare with AIC and BIC.

## Method

- Gaussian likelihood with known sigma; on the same data the additive constant is identical across
  candidates, so only `2k + RSS/sigma^2` (AIC) and `k ln N + RSS/sigma^2` (BIC) are compared. Lower is
  better; the reported quantity is delta-AIC relative to the best candidate.
- Two scenarios, each testing a Phase 1 prediction. Script:
  `experiments/04_discrimination/discriminate.py` (writes `results/discrimination_summary.txt`).
- Important limitation: each cell is a single noise realization. Large delta-AIC values are robust to
  this; small ones (roughly < 3) are within noise-draw variability and are not quantitative.

## Results

Scenario A, truth = M3, candidates M1 vs M3 (delta-AIC of the loser):

| channel       | winner | delta-AIC |
| ------------- | ------ | --------- |
| mRNA-only     | M1     | M3 +2.0   |
| protein-only  | M3     | M1 +0.6   |
| mRNA+protein  | M3     | M1 +25.0  |

Scenario B, M2 vs M3 at symmetric feedback:

| truth | channel       | winner | delta-AIC of loser |
| ----- | ------------- | ------ | ------------------ |
| M2    | mRNA-only     | M2     | M3 +1017           |
| M2    | protein-only  | M2     | M3 +2.5            |
| M2    | mRNA+protein  | M2     | M3 +1139           |
| M3    | mRNA-only     | M2     | M3 +0.4            |
| M3    | protein-only  | M3     | M2 +0.3            |
| M3    | mRNA+protein  | M3     | M2 +26.3           |

## Hypotheses, with evidence for and against

- mRNA cannot separate M1 from M3. Evidence for: under mRNA-only, M1 wins by exactly +2.0, which is the
  AIC parameter penalty for M3's extra parameter with zero fit improvement (their mRNA dynamics are
  identical). Status: confirmed, and the +2.0 is interpretable, not noise.
- Protein observation reveals M3 (the original phrasing). Evidence against: protein-only barely favors
  M3 (+0.6, a tie). It is observing BOTH channels that exposes M3 (+25.0), because then M1 cannot fit
  mRNA and the repressed protein with one parameter set. Status: revised. The discriminating
  information is in the joint mRNA+protein fit, not protein alone.
- The mRNA channel separates M2 from M3. Evidence for and a caveat: when the truth is M2, mRNA crushes
  M3 (+1017), because only M2 regulates mRNA and M3 cannot reproduce that. But the effect is asymmetric:
  when the truth is M3, M2 can mimic it by setting its own mRNA feedback to near zero, so mRNA-only
  cannot reject M2 (+0.4). Only the joint mRNA+protein fit separates them in that direction (+26.3).
  Status: confirmed with an asymmetry. Discrimination power depends on which model generated the data,
  because M2 is the richer model in the mRNA channel.

## Open items (backward verification needed)

- The small delta-AIC cells (0.3 to 2.5) are single-realization and not quantitative. The next step is
  to repeat each scenario over many noise draws and report the distribution of delta-AIC or a model-
  selection error rate, so "weak discrimination" becomes a number rather than a one-shot verdict.
- Feedback strength was fixed (1/kappa = 50 molecules). Discrimination should depend on it; a sweep over
  feedback strength (and over whether M2 is in its oscillatory regime) is untested.
- AIC and BIC agreed on the winner in every cell here; BIC's larger penalty would matter mainly in the
  borderline cells, exactly where replication is needed.

## Notes

- This experiment uses the same simulator and likelihood as Phases 2 and 3, so it inherits the same
  normalization and known-sigma assumptions; a per-channel noise model is a later refinement.
