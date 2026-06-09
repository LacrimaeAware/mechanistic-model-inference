# Handoff brief

A point-in-time snapshot for another agent (or a later session) picking up this project. For the
fixed working agreement read [`AGENTS.md`](AGENTS.md); for the plan read [`docs/roadmap.md`](docs/roadmap.md).
This file records the current state and the open questions, written to be read cold.

## What the project is

Parameter identifiability and inference for small mechanistic ODE models of gene and protein dynamics.
The question is statistical: given a published model and noisy data, which parameters can be estimated,
how well, and from what measurements. The first target is the 2023 negative-autoregulation paper
(three nested mRNA/protein models). See [`docs/literature.md`](docs/literature.md).

## State as of this brief

- Phase 0 (scaffold), Phase 1 (deterministic and stochastic), and Phase 2 (structural map and practical
  profiles for M2 and M3) are done and pushed to `main`. Phase 3 (recovery + MCMC posterior) and Phase 4
  (AIC/BIC discrimination) have initial experiments that still need backward verification.
- The three target models (M1 simplistic, M2 transcriptional, M3 translational) are implemented in
  [`src/mechanistic_inference/models.py`](src/mechanistic_inference/models.py) as the paper's
  dimensionless system, with a closed-form M1 solution, a steady-state solver, and an oscillation
  discriminant.
- Reproduction script: [`experiments/01_reproduce_autoregulation_models/reproduce.py`](experiments/01_reproduce_autoregulation_models/reproduce.py)
  (writes figures and a summary to `results/`, which is gitignored).
- Stochastic SSA: [`src/mechanistic_inference/stochastic.py`](src/mechanistic_inference/stochastic.py)
  with reproduction script `experiments/01_reproduce_autoregulation_models/reproduce_stochastic.py`.
- Phase 2 identifiability map: `make_dimensional_simulator` in `models.py` feeds the carried-over
  Fisher/profile tooling; analysis and figure in `experiments/02_identifiability/identifiability_map.py`,
  write-up in `experiments/02_identifiability/identifiability.md`.
- Tests: 25 passing (6 baseline identifiability, 8 model, 3 stochastic, 8 Phase 2). Run with
  `python -m pytest tests/ -q` (a venv with `requirements.txt` plus `pytest` is assumed; `conftest.py`
  puts `src/` on the path).

## What is and is not established

What the reproduction supports, with the reason:
- M1 numerical integration matches the closed form (Eq. 8) to about 3e-8, which checks the solver.
- Damped oscillation appears only in M2 and only for strong inhibition; M1 and M3 do not overshoot.
  The numerical overshoot agrees in sign with the analytic discriminant `(theta-1)^2 - 4 h1(p*)`.
- Response and return times order as M1 (slowest) then M2/M3, with M3 slightly slower than M2, which
  matches the direction of the paper's Fig. 5.

What is not established, and why to be cautious:
- The match to the paper's figures is qualitative (same shapes and ordering), not a point-by-point
  comparison against digitized figure data. The reproduction has not been checked against the authors'
  own code.
- Two equation details were judgment calls resolved by reading the rendered PDF: the Hill argument
  `(kappa P)^n` and the multiplicative signal entry `(1 + S_a) F(p)`. They are internally consistent
  and match the stated fold-change, but an independent reading could differ.
- The stochastic (Gillespie) reproduction recovers the protein-CV ordering M2 > M1 > M3 and lands
  within about 0.01-0.02 of the paper's protein CVs (0.16/0.14/0.10). The Hill-in-propensity
  formulation is a modeling choice the paper does not write out, so the exact CVs are not claimed to
  match its implementation; the M2 mRNA CV in particular runs higher than reported (about 0.96 vs 0.89).

## Two structural facts logged for later phases

- M3 leaves mRNA unregulated, so M3 mRNA is identical to M1 mRNA. The mRNA channel cannot separate M1
  from M3; that separation needs the protein channel. (Bears on protein-only vs mRNA+protein schemes.)
- M2 and M3 share the same protein steady state when `kappa_tilde_M = kappa_tilde_P` and `n_M = n_P`.
  This is a discrimination degeneracy in the protein-level observable.

These two are hypotheses about where identifiability and discrimination will be hard (the second is a
discrimination question, Phase 4).

## Phase 2 result so far

Structural identifiability (Fisher rank) for the dimensional parameters: in all three models,
protein-only observation leaves exactly one non-identifiable direction, the transcription/translation
rate product (k_m vs k_p); the regulation strength kappa is identifiable; observing mRNA restores full
rank. A hypothesis that transcriptional feedback would break the product degeneracy was refuted, both
analytically (the scaling k_m c, k_p / c leaves the protein trajectory invariant) and via the Fisher
null direction. Practical profiles on M3 agree (k_m flat under protein-only, bounded with mRNA; kappa_P
bounded but weakly). A methodology audit (`verify_methodology.py`) corroborates the structural result
three independent ways (analytic scaling symmetry, Fisher rank, MLE recovery) and shows a noise sweep
where k_m stays non-identifiable at every noise level while kappa crosses identifiable -> not as noise
grows (protein-only between 1-3%, mRNA+protein between 3-10%). Practical profiles are done for both M2
and M3; under protein-only kappa is only marginally identifiable for both (identifiable in 6/10 of noise
draws for M3 kappa_P, 5/10 for M2 kappa_M), and both are clearly identifiable with mRNA. See
`experiments/02_identifiability/identifiability.md`.

## Phase 3 and Phase 4 (initial, need verification)

- Phase 3 (`experiments/03_inference`): recovery on M3 (both channels recover all parameters;
  protein-only recovers the k_m k_p product but not the split) and MCMC posteriors. M1 MCMC is verified
  converged (steps/tau = 138) and shows the protein-only ridge (corr -0.94) that mRNA closes: a fourth
  route to the degeneracy. M3 MCMC does NOT fully converge for protein-only (the flat ridge is
  diffusion-limited; see inference.md for the approaches tried and why each failed), so it is reported
  preliminary with M1 as the anchor.
- Phase 4 (`experiments/04_discrimination`): AIC/BIC. mRNA cannot separate M1 from M3; the joint
  mRNA+protein fit separates models a single channel cannot; M2 vs M3 is asymmetric (mRNA rejects M3
  when the data come from M2, but M2 can mimic M3). Large delta-AIC robust; small cells are single-
  realization and their replication is blocked on fitting speed (below).
- These experiments reproduce by running the scripts (results write to the gitignored `results/`); no
  unit tests were added for them since the fits and MCMC are slow.

## Tooling note (the fitting bottleneck, resolved)

- The original `fit_mle` (Nelder-Mead, unbounded) was slow and fragile: fitting a wrong model sent it
  into stiff parameter regions where individual `solve_ivp` calls were pathologically slow, with no
  per-fit cap. This made the discrimination replication impractical (one run burned ~50 minutes of CPU).
- Fixed by `fit_least_squares` (`mechanistic_inference.identifiability`): trust-region least squares on
  the residuals, with box bounds (theta0 +/- log_bound) that keep the optimizer out of the slow regions
  and a `max_nfev` cap. It is ~8x faster per fit and reliably bounded; the borderline-cell replication
  now runs in ~30 s instead of not finishing. The discrimination scripts use it. The old `fit_mle` is
  kept for the existing profile-likelihood tests.
- Habit adopted: long-running scripts flush per-iteration progress so a run can be monitored live
  instead of buffering silently.

## Data note

Phase 2 does not need an external dataset. Identifiability here is simulation-based: simulate from the
model at known parameters, add synthetic measurement noise, then test what is recoverable. Real
experimental time-course data would only be needed for a later, optional "fit to real measurements"
direction; the target paper provides no such dataset. The sibling project's DREAM datasets are for a
different problem (network-structure inference) and are not required here.

Source PDFs and their text extractions live in `data/raw/` and are gitignored (copyrighted, not for a
public repo). Do not commit them.

## Coordination

If more than one agent works here at once, avoid editing the same files in parallel. Record conclusions
in the relevant experiment note or in this brief rather than only in chat, so the others can read them.
