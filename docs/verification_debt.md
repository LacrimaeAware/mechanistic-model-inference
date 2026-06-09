# Verification debt and targets

Results here are provisional until independently re-verified (see AGENTS.md). This file lists, per
result, the specific things that could be wrong and should be re-checked by a second method or a fresh
code path, not by re-running the same script. Keep it current as items are checked or new results land.

## The standard (for any agent or person continuing this work)

- Before claiming a result, enumerate its verification targets: the specific places it could be wrong
  (the likelihood form, AIC/BIC bookkeeping, optimizer convergence, model code, finite-difference
  accuracy, the noise model, single-realization versus distribution). Check each, and state which were
  checked and which were not.
- Re-verify by a different route (analytic, a second optimizer or sampler, a different parameterization),
  not by re-running the same code.
- Assume the known failure modes recur: optimizer initialization (use multi-start; where a nested-model
  comparison applies, check the nesting floor delta-AIC >= -2), conflated settings between runs (hold
  grid, fitter, and noise fixed when comparing), and single-realization claims (report rates or
  distributions, not one fit).
- Do not call a result "trustworthy", "settled", "verified", or "definitive" on the basis of self-checks.

## On firm ground (tested or analytic cross-check)

- Phase 1 reproduction (models; deterministic and stochastic): correctness tests, analytic M1 check.
- Phase 2 structural identifiability: Fisher rank, the analytic scaling-symmetry proof, and MLE recovery
  agree (three independent routes).

## Provisional, re-verification owed

- Phase 3 inference: recovery numbers come from one optimizer; M1 MCMC convergence from one sampler; the
  M3 MCMC is known under-converged. Targets: a second optimizer for recovery; an independent posterior
  method (or analytic marginal) for the ridge.
- Phase 4 discrimination: the multi-start rates pass the nesting-floor check and agree between two
  scripts, but the AIC/BIC bookkeeping, the likelihood, and the non-nested M2-vs-M3 fits have not been
  independently re-derived. Targets: recompute AIC by hand for a few cells; a likelihood-ratio or
  Bayes-factor cross-check; confirm the M2/M3 multi-start reaches the global optimum (more starts and a
  different optimizer).
- Phase 3 real-data fit (translation model): the fit and the Fisher rank 5/7 are a first pass. Targets:
  the closed-form solution versus a numerical ODE solve; the Fisher null directions versus an analytic
  derivation; profile likelihood (not only Fisher); the delta/gamma label symmetry; single-cell fits
  versus the mean trace.
- Discriminability map (experiment 06, as built): targets noted in its own write-up - the multi-start
  global optimum at each grid cell, the selection rate versus a likelihood-ratio expectation, and
  sensitivity to the grid, noise, and start choices.
