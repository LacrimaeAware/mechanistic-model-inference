# Parameter identifiability and inference for models of gene autoregulation

Repository: `mechanistic-model-inference`. Readable walkthrough with figures:
<https://lacrimaeaware.github.io/mechanistic-model-inference/>.

## In plain terms

A mechanistic model is a set of equations for how a gene's messenger RNA (mRNA) and protein change over
time, including feedback where the protein slows down its own production. A published paper (Ryzowicz &
Yildirim 2023) builds three such models and compares their behavior: (i) no feedback, (ii) the protein
represses its own *transcription* (mRNA production), (iii) the protein represses its own *translation*
(protein production from mRNA).

This project asks the statistical question the paper does not: **given only noisy measurements, which
numbers in the equations can you actually recover, and can you tell the two feedback mechanisms apart?**

The short answer, in three parts:

- **From protein measurements alone you cannot separate the transcription rate from the translation
  rate** — only their product is recoverable. This is structural: no amount of protein data, and no
  fancier fitting method, fixes it, because the information simply is not there.
- **Measuring mRNA fixes it** — and a single well-timed mRNA measurement is enough; dense protein
  sampling cannot substitute.
- **The two feedback mechanisms are indistinguishable from protein**, and only *asymmetrically*
  distinguishable from mRNA: transcriptional feedback is easy to detect (only it changes the mRNA), but
  translational feedback is hard to confirm even with mRNA.

Honest scope: the identifiability facts are established prior art — we confirmed them for these specific
models with standard methods and a hand proof. The mechanism-discrimination answer is the one genuinely
open piece. This is a modest, correct analysis plus a real-data parameter recovery, not a high-impact
result, and it does not claim to be.

## Findings

1. **Protein-only data cannot separate the transcription and translation rates** (only their product
   k_m·k_p is identifiable), in all three models. Fisher rank is 3/4 (no feedback) or 4/5 (with
   feedback); the single non-identifiable direction is exactly the rate product. Established by an
   analytic scaling-symmetry proof, the Fisher rank, and maximum-likelihood recovery (the product is
   recovered to 0.03 in log units; the individual split is not). The regulation strength is identifiable
   in principle, but only marginally from protein at realistic noise (recovered in 5–6 of 10 noise
   draws). [experiment 02; known prior art, confirmed for these models]

2. **Measuring mRNA restores full identifiability, and one informative timepoint suffices.** The
   protein-only deficiency is structural, so no number of protein timepoints removes it; a single mRNA
   measurement closes it, after which every parameter has a finite, modest standard error (~0.02–0.16 in
   log units). Deterministic Fisher analysis. [experiment 07]

3. **The two regulatory mechanisms are indistinguishable from protein; with mRNA, distinguishability is
   asymmetric.** Measured deterministically (how well the wrong mechanism mimics the right one's
   noise-free output): from protein, the best mimic leaves a residual of only 0.02–0.06% of the signal,
   far below realistic noise (1–10%), at every feedback strength. With mRNA, transcriptional
   autoregulation is clearly detectable (the translational model cannot reproduce its mRNA — 6.7–8.5%
   residual), but translational autoregulation is only marginally distinguishable (the transcriptional
   model mimics it to 0.6–2.1%, near the noise floor), easing as feedback strengthens. This is the open
   question the source paper raised but never analyzed. [experiment 09]

4. **On real single-cell data, the pipeline recovers the published degeneracy and real parameters.**
   Fitting the translation model to the Fröhlich et al. (2018) single-cell GFP dataset, the Fisher rank
   recovers the published non-identifiability (the amplitude product; rank 5/7), and the kinetic
   parameters of 200 real cells are recovered (a fast and a slow timescale, ~3.8 h and ~30 h half-lives,
   ~3-fold cell-to-cell variation in expression). A first pass, not independently re-verified.
   [experiment 05; provisional]

## What is solid, provisional, and retracted

- **Solid** (analytic / deterministic / tested): the reproduction of the paper's models (experiment 01),
  the identifiability core (experiment 02, backed by a proof, not just the numerical method), the
  experimental-design result (experiment 07), and the deterministic mechanism distinguishability
  (experiment 09). 33 tests pass, including four that **independently re-derive** the load-bearing
  findings (a symbolic proof of the rate degeneracy among them); see
  [`docs/verification_audit.md`](docs/verification_audit.md).
- **Provisional** (a first pass, re-verification owed): the real-data fits (experiment 05); the
  experimental-design and distinguishability numbers depend on the multi-start fit reaching the global
  optimum. See [`docs/verification_debt.md`](docs/verification_debt.md).
- **Retracted**: the Monte-Carlo "discrimination selection rates" (experiments 04, 06). They compared
  two models that fit protein data near-identically, so the AIC winner was decided by sub-noise wobble
  and swung between runs. Replaced by the deterministic structural statement (finding 3). This is the
  project's main methodological lesson.

## Where this sits in the literature

A verified literature review ([`docs/literature_research.md`](docs/literature_research.md),
[`docs/literature_references.md`](docs/literature_references.md)) places the work honestly:

- The protein-only rate-product non-identifiability and its resolution by observing mRNA are **textbook
  prior art** for this model class (Pieschner, Hasenauer & Fuchs 2022; Fröhlich 2018; Villaverde 2019).
  Our methods (profile likelihood, Fisher information, MCMC) are the standard ones (Raue 2009).
- The source paper (Ryzowicz & Yildirim 2023) builds the M1/M2/M3 models but performs **no
  identifiability analysis, no model discrimination, and no data fitting** — so the mechanism-
  discrimination question (finding 3) is the genuinely open, non-tautological piece.

## What did not pan out (and is not re-attempted)

- **Model discrimination via AIC selection rates** (experiments 04, 06): unstable on near-tied models;
  retracted in favor of the deterministic distinguishability.
- **MCMC for the regulated model M3** (experiment 03): the flat non-identifiable ridge is diffusion-
  limited; the sampler did not converge. The simpler M1 posterior (converged) is the demonstration.
- **Real Hes1 autoregulation data** (experiment 08): the traces are noisy, decay from a high start, and
  need a delayed model with a stochastic noise treatment; the project's non-delay models do not fit them,
  so no fit was forced. Recorded as an honest stop.

## Lessons (methodological, transferable)

1. Do not apply a noisy statistic to a near-degenerate problem and read its fluctuations as results. When
   a number swings between runs, suspect the measurement, not the sampling.
2. Identifiability is a property of the data and model, not of the method: if a parameter is
   non-identifiable, no optimizer, neural network, or Bayesian method can recover it.
3. Prefer deterministic / structural checks over Monte-Carlo rates where possible; they cannot flip-flop.
4. A result is provisional until re-derived by a different route. Passing one's own checks is not
   certification. (See [`AGENTS.md`](AGENTS.md), [`docs/verification_debt.md`](docs/verification_debt.md).)

## Experiment log

| #  | experiment | method | result |
|----|---|---|---|
| 01 | reproduce the models | ODE/SSA simulation vs the paper | oscillation only in the transcriptional model; faster response under feedback; stochastic noise M2 > M1 > M3 (matches the paper); M1 matches its closed form to 3e-8 |
| 02 | identifiability | Fisher rank + profile likelihood + analytic proof | protein-only cannot separate k_m, k_p (rank 3/4, 4/5; null = product); mRNA restores full rank; regulation strength marginal from protein (5–6/10 draws) |
| 03 | inference | MLE recovery + MCMC | protein-only recovers the product not the split; M1 posterior is a ridge (corr −0.93) that mRNA closes; M3 MCMC under-converged |
| 04 | discrimination (AIC) | AIC selection rate | RETRACTED: unstable on near-tied models |
| 05 | real data | fit translation model to Fröhlich 2018 | recovers the published amplitude degeneracy (rank 5/7); 200 single cells fit (half-lives ~3.8 h / ~30 h); provisional |
| 06 | discriminability map (AIC) | AIC rate over feedback × channel | RETRACTED: protein-only rates non-reproducible (sub-noise) |
| 07 | experimental design | deterministic Fisher | protein-only structurally rank-deficient for any protein sampling; one informative mRNA timepoint restores full rank |
| 08 | Hes1 real autoregulation | characterize real Hes1 data | honest stop: noisy, needs a delayed + stochastic model the project does not have |
| 09 | structural distinguishability | deterministic mutual-fit residual vs noise | mechanisms indistinguishable from protein (0.02–0.06%); with mRNA, transcriptional clear (6.7–8.5%), translational marginal (0.6–2.1%) |

## Reproduce

```bash
# Python 3.13, dependencies in requirements.txt; a venv with pytest assumed
.venv\Scripts\python.exe -m pytest tests/ -q                                  # 33 tests
.venv\Scripts\python.exe experiments/02_identifiability/identifiability_map.py # any experiment
.venv\Scripts\python.exe experiments/09_structural_distinguishability/structural_distinguishability.py
```

Datasets (`data/`) and generated figures/tables (`results/`) are git-ignored; the real datasets
(Fröhlich transfection on Zenodo, Hes1 from the GPosc repo) are downloaded on demand by their loaders.
The committed figures used by the website live in `docs/figures/`.

## Layout

```text
mechanistic-model-inference/
├── src/mechanistic_inference/   # models and the identifiability/inference tooling
├── experiments/                 # one folder per experiment, each with a write-up and a script
├── docs/                        # summary, methodology review, literature, the website (docs/index.html)
└── tests/                       # 33 correctness tests
```

## Further reading

- [`docs/summary.md`](docs/summary.md): the one cohesive story across all experiments.
- [`docs/methodology_review.md`](docs/methodology_review.md): whether each method was the right tool.
- [`docs/literature_research.md`](docs/literature_research.md) and
  [`docs/literature_references.md`](docs/literature_references.md): the verified literature and novelty
  verdict, with citations.
- [`docs/verification_audit.md`](docs/verification_audit.md): the independent re-derivation of the load-bearing findings and the per-experiment methodology audit.
- [`docs/verification_debt.md`](docs/verification_debt.md): what is verified and what still owes a check.
- Each `experiments/NN_*/` directory has the full numbers and the per-experiment write-up.
