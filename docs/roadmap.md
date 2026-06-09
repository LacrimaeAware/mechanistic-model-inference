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

## Phase 3: inference under noise (initial experiments done; verification pending)

- Maximum-likelihood fitting on simulated noisy data; recovery as a function of noise and sampling.
  Done for M3 (`experiments/03_inference`): both channels recover all parameters; protein-only recovers
  the k_m k_p product but not the split.
- Then Bayesian / MCMC (emcee) for posterior uncertainty. Done for M1: the protein-only posterior is a
  ridge along log k_m + log k_p = const (correlation -0.93); observing mRNA closes it into a bounded
  blob. MCMC on the regulated models and formal convergence diagnostics are deferred.
- Deliverable: recovery and uncertainty results, honest about what cannot be recovered. Initial version
  done; a sampling-density sweep and M3 MCMC remain.

## Phase 4: model discrimination and experimental design (initial experiments done; verification pending)

- Whether the autoregulation modes can be told apart from data (AIC/BIC), and under what measurement
  scheme. Initial experiments (`experiments/04_discrimination`): mRNA cannot separate M1 from M3; the
  joint mRNA+protein fit separates models where a single channel cannot; M2 vs M3 discrimination is
  asymmetric (mRNA rejects M3 when the truth is M2, but M2 can mimic M3). Large delta-AIC results are
  robust; small ones need replication over noise draws.
- Optimal experimental design: which added measurement most improves identifiability or discrimination.
  The results already point at observing mRNA as the high-value measurement; a systematic design metric
  is not yet computed.

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
- Phase 2 finished (M2 and M3 practical profiles done). Phase 3 and Phase 4 have initial experiments:
  recovery and an MCMC posterior ridge (Phase 3), and AIC/BIC discrimination (Phase 4). These are
  exploratory and partly need backward verification (noise-draw replication for the borderline
  discrimination cells; M3 MCMC; a feedback-strength sweep).
- Fitting bottleneck fixed: a bounded least-squares fitter (`fit_least_squares`) made large sweeps
  feasible. The borderline discrimination cells were replicated over noise draws (28 tests passing);
  the replication corrected a single-run overstatement (protein-only distinguishes M1 from M3 in 88% of
  draws, not the one-shot tie reported earlier). The M3 protein-only MCMC is diffusion-limited and
  reported preliminary; the converged M1 MCMC is the quantitative ridge anchor.
- Real data acquired: the Frohlich et al. 2018 single-cell GFP translation dataset (Zenodo), loaded and
  visualized (`experiments/05_real_data`). It is the protein-observed / mRNA-latent identifiability
  question on real data, with a published benchmark (Pieschner et al. 2022) to check against.
- Next: implement the translation model and fit it to the real data, then run identifiability on the
  real fit and compare to the published degeneracy. Stretch: a discriminability map (feedback strength
  x noise) now that fitting is fast.
