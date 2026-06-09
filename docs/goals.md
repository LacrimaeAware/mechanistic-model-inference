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
