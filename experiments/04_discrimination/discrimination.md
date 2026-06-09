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

## Replication over noise draws (backward verification, done) and a fitting bug it exposed

The borderline cells were replicated over 25 noise draws each
(`experiments/04_discrimination/replicate_discrimination.py`). A controlled, fully fixed-settings check
of the most contested cell (`m1_m3_controlled.py`) exposed a real bug first: the model fit is
**initialization-dependent**, because near kappa = 0 the protein trace barely depends on kappa (a flat
gradient), so a single-start optimizer gets stuck wherever it began. A diagnostic that exploits the
nesting of M1 inside M3 (so M3 can never fit worse, delta-AIC >= -2) showed this directly: a single
start violated that floor in 5 of 60 draws (one start) or 42 of 60 (a different start). The fix is
**multi-start over kappa** (several kappa values plus the simpler model's fit), taking the global best;
that brings floor violations to ~0 and makes the AIC comparison reflect the real optimum.

Multi-start results (the trustworthy ones; "select rate" = fraction of draws the true model wins on AIC):

| comparison              | channel       | true-model select rate | verdict            |
| ----------------------- | ------------- | ---------------------- | ------------------ |
| M1 vs M3 (truth M3)     | mRNA-only     | 0%                     | indistinguishable (parsimony picks M1) |
| M1 vs M3 (truth M3)     | protein-only  | 88% (controlled: 85%)  | mostly distinguishable |
| M2 vs M3 (truth M2)     | protein-only  | 88%                    | mostly distinguishable |
| M2 vs M3 (truth M3)     | mRNA-only     | 24%                    | indistinguishable (M2 mimics M3) |
| M2 vs M3 (truth M3)     | protein-only  | 44%                    | indistinguishable (coin flip)    |

What the bug fix changed (the single-start numbers were contaminated, so treat the earlier table and the
single-run `discriminate.py` output as superseded):
- M2 vs M3, truth M2, protein-only: 96% -> 88% (the wrong model M3 can actually fit better than the stuck
  optimizer found, so M2 is favoured less often).
- M2 vs M3, truth M3, protein-only: 68% "marginal" -> 44% "indistinguishable". This is a changed
  conclusion, not a tweak: from protein alone you essentially cannot tell M2 from M3 when the data come
  from M3.

Reliable headline (asymmetric discrimination):
- M1 vs M3: mRNA alone cannot separate them (identical mRNA); protein alone mostly can (~85-88%).
- M2 vs M3: when the truth is M2 you can mostly tell (protein 88%, mRNA ~100% since only M2 regulates
  mRNA); when the truth is M3, M2 can mimic it, so neither single channel is reliable (mRNA 24%,
  protein 44%).

## Open items

- Discrimination depends on feedback strength, held fixed here at 1/kappa = 50 molecules; a sweep over
  feedback strength (and over M2's oscillatory regime) is the natural next step (the "discriminability
  map"), now feasible with the fast multi-start fitter.
- BIC (heavier penalty than AIC) would shift the borderline rates; worth reporting alongside the sweep.
- A per-channel noise model (mRNA and protein measured at different precision) is a later refinement.
- The single-realization `discriminate.py` table is kept for the record but is superseded by the
  multi-start replication; single-start AIC on these models is not reliable.

## Notes

- This experiment uses the same simulator and likelihood as Phases 2 and 3, so it inherits the same
  normalization and known-sigma assumptions; a per-channel noise model is a later refinement.
