# Research questions

Grounded in the first target (the 2023 negative-autoregulation models; see docs/literature.md) and the
inference toolkit. These will grow and sharpen as the work proceeds.

## Reproduction

- Can the three published ODE models (no regulation; transcriptional negative autoregulation;
  translational negative autoregulation) be reproduced from their equations and parameters, matching
  the paper's deterministic behavior (faster response and quicker return under regulation; oscillation
  only in the transcriptional model)?
- Can the paper's stochastic-noise ordering (transcriptional noisiest, simplistic next, translational
  least) be reproduced with a standard stochastic simulation (for example Gillespie)?

## Identifiability (the core question)

- For each model, which parameters are structurally identifiable (Fisher-information rank) from
  (a) protein-only and (b) mRNA-plus-protein time-course observation?
- Which parameters are practically identifiable (bounded profile likelihood) at realistic noise and
  sampling density?
- The no-regulation baseline is known: transcription and translation rates are non-identifiable from
  protein alone and identifiable when mRNA is observed. Do the regulation parameters introduce new
  non-identifiable combinations?

## Model discrimination

- Can transcriptional and translational negative autoregulation be distinguished from time-course data
  (likelihood ratio, AIC/BIC), and under what measurement scheme and noise level?
- Is there a regime where the two mechanisms are statistically indistinguishable despite producing
  different dynamics?

## Experimental design

- Which single added measurement (observing mRNA, finer time resolution, a changed input step) most
  improves identifiability or model discrimination?

## Orthogonal and transfer questions (open, to be scoped honestly)

- The inference-for-dynamical-systems toolkit (likelihood, Fisher information, profile likelihood,
  filtering) is domain-general. Which parts have established applications outside molecular biology
  (finance state-space and filtering models, image or other time-series), and which are only plausible
  bridges? This must be checked against the literature before any transfer is claimed.
