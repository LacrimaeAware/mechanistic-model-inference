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

- Phase 0 (scaffold), Phase 1 (deterministic and stochastic), and the Phase 2 structural
  identifiability map are done and pushed to `main`. The Phase 2 practical sweep is partial.
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
- Tests: 24 passing (6 baseline identifiability, 8 model, 3 stochastic, 7 Phase 2). Run with
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
bounded but weakly). The remaining Phase 2 work is a practical sweep over noise and sampling, and
profiles for M2. See `experiments/02_identifiability/identifiability.md`.

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
