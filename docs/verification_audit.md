# Verification audit

**Scope:** all 9 experiments and the shared library (`src/mechanistic_inference/`). An independent
end-to-end review on two axes: **program correctness** (does the code compute what it claims) and
**design validity** (does the experiment actually answer its stated question, and is the method
research-backed and free of flaws). The load-bearing findings were re-derived by a *different route* than
the one that produced them, so the conclusions do not rest on a single implementation.

## Independent re-derivations of the load-bearing findings

These are committed as `tests/test_independent_verification.py` (4 tests, all pass) so any reader can
re-run them. Each uses a method independent of the original code.

1. **Protein-only rate degeneracy (Findings 1 and 2), proven symbolically.** Using `sympy` and nothing
   from this project's tooling, the mRNA->protein cascade was solved in closed form and differentiated:
   the protein's sensitivity to `log k_m` and to `log k_p` are *identical* (their difference simplifies
   to exactly 0), and each equals the protein itself, so the protein is proportional to the product
   `k_m·k_p` and the two rates have exactly one non-identifiable direction. The mRNA's sensitivity to
   `log k_p` is *exactly 0* while its sensitivity to `log k_m` is nonzero, so observing mRNA carries
   information about `k_m` and breaks the tie. This is an algebraic proof, not a numerical fit.
2. **The degeneracy holds for the regulated nonlinear models.** A high-precision integration (tolerance
   1e-11) confirms that scaling `k_m -> c·k_m, k_p -> k_p/c` leaves the M2 and M3 protein trajectories
   invariant to ~1e-11 (machine precision), so the degeneracy is not an artifact of the linear M1.
3. **The Fisher null direction equals the analytic one** (`|null · [1,0,-1,0]/sqrt2| = 1.000000`), so the
   numerical Fisher rank agrees with the proof.
4. **The distinguishability result is the global optimum, not an under-fit.** Doubling the optimizer
   starts (5 -> 11) does not change the best mimic residual (0.0617% either way), so the small residual
   in Finding 3 is real, not a failure to fit.

Conclusion of the re-derivation: Findings 1, 2, and 3 are confirmed by independent methods; Finding 1 in
particular now rests on a symbolic proof.

## Per-experiment audit

**01 reproduce the models.** Correctness: the dimensionless models match the paper's equations
(confirmed earlier against the rendered PDF), and M1 matches its closed-form solution to 3e-8; tests
cover the models and the stochastic ordering. Design: forward simulation compared to the paper's reported
behavior is unambiguous, and the reproduced behavior (oscillation only in the transcriptional model,
faster response under feedback, noise M2 > M1 > M3) matches. Flaw/caveat: the comparison to the paper's
figures is qualitative (shapes and orderings, not digitized values), and the stochastic propensity for
the Hill term is a modeling choice the paper does not write out. Status: **sound**.

**02 identifiability.** Correctness: Fisher rank, the null eigenvector, and profile likelihood are
computed correctly; the spectral gap is clear (the null eigenvalue is ~1e-16 of the largest, the next is
~1e-5), so the rank is not a threshold judgment call. Design: Fisher rank is the standard structural-
identifiability proxy, and here it is backed by the *symbolic proof* above, so the conclusion is method-
independent and global, not merely a local numerical result. Research-backed: the product non-
identifiability and its resolution by observing mRNA are established prior art (Pieschner 2022; Villaverde
2019), and profile likelihood is the canonical practical-identifiability tool (Raue 2009). Flaw/caveat:
the *practical* identifiability of the regulation strength from protein is realization-dependent
(identifiable in 5-6 of 10 noise draws), correctly reported as marginal. Status: **sound** (the strongest
part).

**03 inference.** Correctness: MLE recovery and the MCMC are computed correctly; the M1 posterior ridge
(correlation -0.93) reproduces the degeneracy. Design: the recovery agrees with experiment 02, a
consistency check across methods. Flaw/caveat: the regulated-model (M3) MCMC does not converge on the
flat ridge (steps / autocorrelation time below threshold) and is reported as preliminary, not used for a
claim. Status: **sound for M1; the M3 posterior is an acknowledged gap**, consistent with the literature
(ensemble MCMC is poor on flat ridges).

**04 and 06 discrimination via AIC selection rate — RETRACTED.** Correctness: the AIC arithmetic is
correct. Design: **invalid for the question.** The two models fit protein-only data to under 1% of peak
(experiment 09), so an AIC comparison is decided by sub-noise wobble; the resulting "selection rate" was
not reproducible (two of our own runs gave 16% and 80% for the same cell). This is the project's main
design-validity failure: a noisy statistic applied to a near-degenerate problem. It was caught (the rate
swinging between runs), diagnosed (the deterministic mutual-fit check), and replaced by experiment 09.
Status: **retracted; superseded by experiment 09.**

**05 real data.** Correctness: the closed-form translation model fits the real Frohlich mean trace to
1.4% of peak; the Fisher rank of the separated 7-parameter model is 5/7 with the two null directions in
the scale/k/m0 block. Design: the structural-identifiability finding recovers the *published* result
(Pieschner et al. 2022, the amplitude is a product), which is external corroboration that the method is
applied correctly on real data. The single-cell recovery fits 200 cells with sensible kinetics. Flaw/
caveat: it fits the population *mean* (not a single-cell mixed-effects model), the two decay-rate labels
are not identifiable (a known symmetry), only Fisher (not profile) identifiability was run, and the
parameter estimates have no per-cell confidence intervals. Status: **provisional, externally
corroborated.**

**07 experimental design.** Correctness: the deterministic Fisher rank and the standard errors are
computed correctly; the minimum-mRNA-timepoints search is correct after a fix (an early version reported
2 because the evenly-spaced first point landed at t=0 where mRNA is identically zero; placed in the
informative window it correctly reports 1). Design: deterministic, so stable; the claim that protein-only
is rank-deficient for *any* protein sampling is exactly the established principle that structural
unidentifiability cannot be removed by more data of the same observable (Villaverde 2019). Flaw/caveat:
the Fisher matrix is finite-difference (accuracy), and "one timepoint restores full rank" is a
*structural* statement; the practical standard error with a single mRNA point is not quantified. Status:
**sound**.

**08 Hes1 real autoregulation data.** Correctness: the data is loaded with correct per-cell NaN handling
(an early version mis-averaged over NaNs; fixed). Design: this is an **honest stop, not a result.** The
project's non-delay models do not match real Hes1 (the traces decay from a high start while the models
rise; sustained oscillation needs a transcriptional delay; the noise needs a stochastic treatment, which
the source paper handled with Gaussian processes). No fit was forced. Status: **recorded honestly as a
non-result.**

**09 structural distinguishability.** Correctness: the deterministic mutual-fit residuals are computed
correctly, and the global optimum is confirmed (re-derivation check 4: more starts do not improve the
fit). Design: deterministic, so it cannot flip-flop; it is the correct replacement for the retracted AIC
rates. The asymmetry (transcriptional detectable from mRNA, translational only marginally) is consistent
with the AIC-era qualitative hints and with the literature's note that the autoregulated case is subtler.
Research-backed: this is the open question the source paper raised but never analyzed. Flaw/caveat: the
residual uses per-channel normalization (an absolute per-channel noise model would rescale the comparison
to the noise lines). Status: **sound** (the contribution).

## Cross-check of reported numbers

Every quantitative claim in the README, the website (`docs/index.html`), and `docs/summary.md` was
written from the saved result files in `results/` (the `*_summary.txt` files) and matches them. The
deterministic experiments (01, 02, 07, 09) were re-run during this audit and regenerate their numbers and
figures cleanly; the full suite is 33 tests passing.

## What could not be independently re-derived

- The **stochastic noise ordering** (experiment 01) and the **MCMC** (experiment 03) are Monte-Carlo; they
  were re-run rather than re-derived, and the M3 MCMC is openly under-converged.
- The **real-data fits** (experiments 05, 08) depend on the external datasets, which are downloaded on
  demand by their loaders rather than committed; the structural-identifiability claim on real data is
  corroborated by the published Pieschner result, but the single-cell parameter estimates have not had an
  independent re-derivation.

## Conclusion

The shared library is correct and the identifiability core is now backed by an independent symbolic proof.
Of the 9 experiments: 5 are sound (01, 02, 07, 09, and 03 for M1), 1 is provisional but externally
corroborated (05), 1 is an honest non-result (08), the M3 MCMC is an acknowledged gap (03), and 2 are
retracted for a design flaw that was caught and corrected (04, 06). No fabricated numbers, no metric
misuse beyond the retracted AIC-rate (which is documented as such), and no circular reasoning were found.
The main lesson the audit confirms: the one genuine methodological error (AIC discrimination on near-tied
models) was a *design*-validity failure, not a coding bug, and it was the deterministic re-derivation,
not more re-running, that exposed and fixed it.
