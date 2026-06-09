# Goals

## Learning goals

The project is also a way to learn the inference toolkit on a real model, in order:

- ODE simulation: integrate small systems, read off steady states, understand solver tolerances.
- Parameter meaning: what each constant controls physically (production rate, degradation rate,
  regulation strength, Hill exponent, input gain).
- Least-squares and maximum-likelihood fitting: fit a model to noisy time-course data, and the
  relationship between the two under Gaussian noise.
- Fisher information: local curvature of the log-likelihood; a near-zero eigenvalue is a
  non-identifiable parameter direction.
- Profile likelihood: practical identifiability and confidence intervals from re-optimized profiles.
- Later: Bayesian / MCMC inference for posterior uncertainty, and particle filtering for stochastic
  or partially observed state-space models.

## Project goals

- Reproduce a published mechanistic model from its equations and match its reported behavior.
- Report which parameters are structurally and practically identifiable, under which measurement
  scheme, with the non-identifiable directions named explicitly.
- Attempt a model-discrimination result: whether transcriptional and translational negative
  autoregulation can be distinguished from data, and what measurement separates them.
- Optimal experimental design: what observation renders an unidentifiable parameter identifiable, or
  discriminates the models.

## What a minimum positive result looks like

A correct, regeneratable identifiability map for one published model: which parameters are estimable
from which data, with the unidentifiable combinations stated. This is a recognized contribution type
in mathematical biology and does not require beating a benchmark. The bar is correct and useful,
verified against a known case (the no-regulation mRNA-to-protein baseline) before any real claim.

## Realistic versus stretch

- Realistic (one person, months): reproduce one model, full structural and practical identifiability,
  maximum-likelihood inference under noise.
- Stretch: Bayesian / MCMC with uncertainty quantification; model discrimination plus optimal
  experimental design; a short write-up.
- Out of scope unless it earns a place: new wet-lab data; large model families; anything that depends
  on beating a benchmark.

## Positioning and relevance (honest)

Stated plainly so the work is aimed correctly:

- This is a small analysis-and-learning project, sibling to stable-grn-inference. It currently sits
  between a learning exercise and a modest analysis contribution. Most results so far are on synthetic
  data and confirm expected structure; the non-obvious finding (the asymmetric mechanism
  discriminability) is the candidate seed of something more.
- The realistic research output is a short, correct note: an identifiability-and-discriminability
  analysis of the target's negative-autoregulation models - which parameters are estimable from which
  measurements, and when the two regulatory mechanisms can be told apart from data and what to measure.
  A recognized but modest contribution type; it does not depend on beating a benchmark.
- Relevance to the target's research: the discriminability/identifiability of the target's specific
  mechanisms is the complementary question the source paper raises (it compares transcriptional versus
  translational autoregulation) but does not answer (whether data can distinguish them). That is the
  relevant core. The real-data work on the open-loop translation dataset validates the tooling but is
  tangential to the target's autoregulation biology; connecting to a real autoregulation system (for
  example Hes1) would make it directly relevant.
- What would make it matter, short of effort alone: a clean, correct, new result about the target's
  models that the source did not produce - chiefly the mechanism-discrimination map plus an
  experimental-design recommendation - reported with explicit limits. A reproduction, or a confirmation
  of a known degeneracy, is competence, not contribution.
- Honest gap: the project studies the statistical properties (identifiability) of the models, which is
  adjacent to but not the same as the dynamical and biological story the source cares about; engaging
  the dynamics more would deepen the understanding of the source work itself.
