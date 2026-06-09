# Literature references (annotated)

The key papers from the verified literature research, organized by topic, so the field is available
without re-researching. Confidence tags reflect the adversarial verification (most load-bearing claims
were confirmed 3-0). The synthesis and the novelty verdict are in `docs/literature_research.md`; the
biological-model sources for the target are in `docs/literature.md`.

## The target paper

- Ryzowicz & Yildirim (2023), "Differential roles of transcriptional and translational negative
  autoregulations in protein dynamics," Molecular Omics 19(1):60-71.
  https://pubs.rsc.org/en/content/articlehtml/2023/mo/d2mo00222a
  Builds exactly the project's M1/M2/M3 two-equation mRNA-protein Hill model. Studies only the dynamical
  roles (stability, oscillation, noise/CV, response/return time) with literature-fixed E. coli parameters.
  Performs NO identifiability analysis, NO statistical model discrimination, NO real-data fitting
  (verified 3-0). This absence is where the project's only open question lives.

## Identifiability of mRNA -> protein models (the "product degeneracy")

- Pieschner, Hasenauer & Fuchs (2022), "Identifiability analysis for models of the translation kinetics
  after mRNA transfection," J Math Biol 84:56. https://link.springer.com/article/10.1007/s00285-022-01739-x
  (PMC9110294). Shows for the translation-after-transfection ODE that the translation rate, initial mRNA,
  and fluorescence scale are identifiable only as a PRODUCT, attributing the degeneracy to Frohlich 2018.
  Also fits the model to real HuH7 transfection data and finds the SDE gives better identifiability than
  the ODE. This is the published result experiment 05 reproduces. [high]
- Frohlich, Reiser, Theis, Hasenauer et al. (2018), "Multi-experiment nonlinear mixed effect modeling of
  single-cell translation kinetics after transfection," npj Syst Biol Appl 4:42.
  https://www.nature.com/articles/s41540-018-0079-7. Origin of the product degeneracy for this model and
  the single-cell GFP dataset used in experiment 05. [high]
- Burton, Manning, Papalopulu, Kursawe et al. (2021), "Inferring kinetic parameters of oscillatory gene
  regulation from single cell time-series data," J R Soc Interface 18:20210393.
  https://royalsocietypublishing.org/rsif/article/18/182/20210393 (PMC8479358). Fits a stochastic Hes5
  auto-NEGATIVE-feedback model to REAL Venus::HES5 single-cell data using a delay-adapted Kalman filter in
  a MALA sampler. Finds the transcription rate poorly inferred from protein only (posterior up to 5x off /
  equal to the prior), with uncertainty cut >50-60% once mRNA (smFISH-type) is added. IMPORTANT: this is a
  DELAYED nonlinear feedback oscillator, and it finds translation well-inferred but transcription poorly
  inferred - a related but DISTINCT non-identifiability from the clean linear product symmetry. [high]

## Identifiability methods and theory (standard toolkit)

- Raue, Kreutz, Timmer et al. (2009), "Structural and practical identifiability analysis ... using the
  profile likelihood," Bioinformatics 25(15):1923. https://academic.oup.com/bioinformatics/article/25/15/1923/213246
  The canonical profile-likelihood reference (~thousands of citations). Structural non-identifiability =
  flat profile (a constant-chi-squared manifold from redundant parameterization); practical = finite
  minimum with an unbounded confidence interval; adding observables can make parameters identifiable. [high]
- Villaverde (2019), "Observability and structural identifiability of nonlinear biological systems,"
  Complexity 2019:8497093. https://onlinelibrary.wiley.com/doi/10.1155/2019/8497093. Structural
  unidentifiability arises from symmetries and "cannot be removed ... unless the new data involves
  modifying the output" - sampling more densely or longer is "doomed to fail." [high]
- Massonis, Villaverde & Banga (2023), "AutoRepar ... structural identifiability and observability,"
  PLOS Comput Biol 19(3):e1011014. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011014
  Symmetries make a model overparameterized hence structurally unidentifiable; automated reparameterization
  (AutoRepar in STRIKE-GOLDD 4.0) removes scaling symmetries. [high; note AutoRepar is niche]
- Villaverde, Pathirana, Frohlich, Hasenauer, Banga (2022), "A protocol for dynamic model calibration,"
  Briefings in Bioinformatics. https://pmc.ncbi.nlm.nih.gov/articles/PMC8769694/. Step-by-step protocol;
  remedy for non-identifiability is new observables or experimental conditions. [high]
- Gutenkunst, Sethna et al. (2007), "Universally sloppy parameter sensitivities in systems biology
  models," PLOS Comput Biol 3(10):e189. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2000971/. Hessian
  eigenvalue spectra span >10^6 across 17 models ("sloppiness"). CAVEAT: sloppiness is NOT
  non-identifiability - see Apgar et al. (2010) "Sloppy models can be identifiable" and Chis, Villaverde,
  Banga, Balsa-Canto (2016) "On the relationship between sloppiness and identifiability." Assess
  identifiability with structural/profile tools, not the eigenvalue spectrum. [high]

## Tools for rigorous (global) structural identifiability

- STRIKE-GOLDD (MATLAB), https://github.com/afvillaverde/strike-goldd ; SIAN (Julia/Maple); 
  StructuralIdentifiability.jl (Julia); DAISY; COMBOS; GenSSI. These give gold-standard global structural
  identifiability via differential algebra or Lie symmetries - stronger than the project's finite-
  difference Fisher rank, and the recommended way to confirm the M2/M3 result.

## Inference on real single-gene data (prior art for "fit real data")

- Burton et al. (2021), above - real Venus::HES5 single-cell time series.
- Pieschner et al. (2022), above - real HuH7 transfection GFP data.
- Heron & Finkenstadt (2007), "Bayesian inference for dynamic transcriptional regulation; the Hes1
  system," Bioinformatics 23(19):2596. https://academic.oup.com/bioinformatics/article/23/19/2596/185373
  Early Bayesian inference on Hes1 dynamics. [medium]
- Calderazzo et al. (2019) - Kalman-filter inference for gene expression, extended by Burton 2021.

## The open question / gap (the only non-tautological route)

No surveyed paper performs statistical model discrimination between transcriptional and translational
negative autoregulation (M2 vs M3) from time-course data, nor fits such a NAR model to a real dataset with
BOTH mRNA and protein co-measured over the same time course. Candidate co-measured systems named:
JAK-STAT CIS/SOCS (Bachmann et al. 2011), Hes/Her oscillators. (Absence in this search is not proof of
absence; a targeted search is warranted before claiming the gap is fully open.)

## Do NOT rely on these (claims that FAILED verification in the research)

- "Structural identifiability is a special case of observability via the augmented-state rank condition"
  as stated (refuted 0-3 in this run; the general idea exists in Villaverde's work but the specific claim
  did not verify).
- A specific "k1*k2*x0 textbook product example" and a "permutation symmetry (k3,k4) example" from one
  source (1-2, not confirmed).
- A "single-species stochastic snapshot" identifiability analogue (0-3; it is a different phenomenon).
