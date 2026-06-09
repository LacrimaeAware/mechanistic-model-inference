# Roadmap

A phased plan from the first concrete target outward. It is expected to evolve as the sibling project
did: start with one model, widen as ideas are tested and pruned. Realistic versus stretch scope is in
docs/goals.md.

## Phase 0: scaffold (done)

- Project structure, conventions (AGENTS.md), goals, research questions, literature.
- Carried-over identifiability tooling, validated on the textbook mRNA-to-protein model.

## Phase 1: reproduce the target models (deterministic part done)

- Transcribe the three ODE models from the paper (no regulation; transcriptional negative
  autoregulation; translational negative autoregulation), exact right-hand sides and parameters, into
  `mechanistic_inference.models`. Done: the dimensionless models (Eqs. 4-7) are implemented and
  covered by `tests/test_models.py`.
- Match the published deterministic behavior: faster response and quicker return under regulation, and
  oscillation in the transcriptional model. Done: reproduced in
  `experiments/01_reproduce_autoregulation_models` (Figs. 4-5 ordering; M1 matches the Eq. 8 closed
  form to 3e-8).
- Optional: reproduce the stochastic-noise ordering with a Gillespie simulation. Done: exact SSA in
  `mechanistic_inference.stochastic` reproduces the protein-CV ordering M2 > M1 > M3
  (`reproduce_stochastic.py`).
- Deliverable: a script and note that regenerate the paper's figures from the equations. Done.

## Phase 2: identifiability (structural map done; practical sweep in progress)

- For each model and each observation scheme (protein-only, mRNA-plus-protein), compute Fisher rank
  and profile likelihoods. Done structurally for all three models and both schemes; practical profiles
  done for M3, partial for the others (`experiments/02_identifiability`).
- Name the non-identifiable directions; compare to the known no-regulation baseline. Done: the only
  protein-only null direction is the k_m vs k_p product, the same in all three models; kappa is
  identifiable; mRNA restores full rank.
- Deliverable: an identifiability map per model. Done (structural); a noise and sampling sweep for the
  practical map is the remaining piece.

## Phase 3: inference under noise

- Maximum-likelihood fitting on simulated noisy data; recovery as a function of noise and sampling.
- Then Bayesian / MCMC (for example emcee) for posterior uncertainty on the identifiable parameters.
- Deliverable: recovery and uncertainty results, honest about what cannot be recovered.

## Phase 4: model discrimination and experimental design

- Whether transcriptional and translational negative autoregulation can be told apart from data
  (likelihood ratio, information criteria), and under what measurement scheme.
- Optimal experimental design: which added measurement most improves identifiability or discrimination.

## Phase 5: extensions

- A second published model (for example a cellular-adaptation model with integral or antithetic
  feedback; see docs/literature.md) once the pipeline is proven on the first.
- Orthogonal transfer of the toolkit to non-biology dynamical data, scoped honestly against the
  literature (see docs/research_questions.md).

## Directions considered (kept, not erased)

- A. Identifiability and inference for small mechanistic models (this roadmap). First focus. Math and
  statistics fit, individual-scale, low confound risk; the bar is a correct and useful analysis, not a
  benchmark win.
- B. Non-additivity and epistasis on combination-perturbation data (a separate, benchmark-heavy
  direction documented in the sibling project). Kept as a backup, not pursued at the same time, to
  keep a one-person effort focused.

## Build status

- Identifiability pipeline built and validated on the textbook mRNA-to-protein model: observing
  protein only, transcription and translation rates are non-identifiable (Fisher rank 3 of 4);
  observing mRNA as well makes them identifiable (rank 4 of 4). See `tests/test_identifiability.py`.
- Phase 1 done (deterministic and stochastic): the three target models are in
  `mechanistic_inference.models`, the exact SSA in `mechanistic_inference.stochastic`, and the paper's
  deterministic behavior and stochastic noise ordering are reproduced
  (`experiments/01_reproduce_autoregulation_models`). Full suite is 17 tests passing (6 identifiability,
  8 models, 3 stochastic). Two structural degeneracies were noted for Phase 2: mRNA cannot separate M1
  from M3, and M2/M3 protein steady states coincide under symmetric feedback.
- Phase 2 structural map done: in all three models the only protein-only non-identifiable direction is
  the k_m vs k_p product; kappa is identifiable; mRNA restores full rank (`experiments/02_identifiability`,
  24 tests passing). A hypothesis that transcriptional feedback breaks the degeneracy was refuted.
- Next: the practical identifiability sweep (noise and sampling) to finish Phase 2, then Phase 3
  (maximum-likelihood recovery, then Bayesian / MCMC).
