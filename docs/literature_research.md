# Literature research: is this novel, and what is the real route forward?

A multi-source, adversarially-verified literature review (21 of 25 claims confirmed 3-0). Bottom line:
the identifiability core the project rests on is established prior art, not novel; the one genuinely open
question is mechanism discrimination, which the source paper raises but never analyzes.

## What is already known (so, not a contribution to re-derive)

- The product-only non-identifiability (in a two-state mRNA->protein model the synthesis/translation rate
  constants are identifiable only as a product from protein-only data, and observing mRNA restores
  identifiability) is established, near-textbook prior art for this model class. Pieschner, Hasenauer &
  Fuchs (2022, J Math Biol 10.1007/s00285-022-01739-x) show it for the translation-after-transfection
  ODE and attribute the product degeneracy to Frohlich et al. (2018, npj Syst Biol Appl). Villaverde
  (2019) and Massonis-Villaverde-Banga (2023, PLOS Comput Biol) give the general scaling-symmetry reason,
  including that such non-identifiability cannot be removed by more/denser/longer data of the same
  observable - only by measuring an extra state (mRNA) or fixing parameters.
- Our methods are the standard, correct ones: profile likelihood (Raue et al. 2009, Bioinformatics, the
  canonical reference), Fisher information, Bayesian/Kalman-MCMC (Burton et al. 2021). Caveat from the
  literature to heed: "sloppiness" (Gutenkunst 2007) is NOT the same as non-identifiability (Apgar 2010;
  Chis et al. 2016) - identifiability must be assessed by structural/profile tools, which we did.
- Real-data fits of small mRNA->protein (auto)regulation models exist: Burton et al. (2021, J R Soc
  Interface) fit a stochastic Hes5 auto-negative-feedback model to real Venus::HES5 single-cell data;
  Pieschner/Fuchs fit unregulated translation to real HuH7 transfection data (the dataset experiment 05
  uses).

So our Phase 2/3/5 results reproduce/confirm known structure with standard methods. Correct and a good
exercise, but not novel.

## The genuinely open question (the only non-tautological route)

Ryzowicz & Yildirim (2023) build exactly the project's M1/M2/M3 models but study only the dynamical roles
(stability, oscillation, noise, response time) with literature-fixed parameters. They perform NO
identifiability analysis, NO statistical model discrimination, and NO real-data fitting. So the open
question is: can transcriptional vs translational negative autoregulation be statistically distinguished
from time-course data, and what measurements are required? The surveyed literature has not settled this
for this mechanism pair.

We have strong preliminary evidence (experiment 04/06, the deterministic structural check) that the
answer is: not from protein alone (the two models fit protein-only data to under 1% RMS, below realistic
noise), and yes once mRNA is observed. Done rigorously, that is a clean, correct answer to an actually-
open question - a modest but genuine contribution.

## Caveats the research flags about our own claims (verification targets)

- The clean linear "only the product" result is proven for UNREGULATED translation. For the AUTOREGULATED
  Hill models (M2/M3) the autoregulated case is "related but distinct" in the literature (Burton finds
  translation well-inferred but transcription poorly inferred under nonlinear delayed feedback). For our
  specific non-delay M2/M3 the scaling symmetry was proven by hand and confirmed by Fisher rank, but a
  rigorous structural-identifiability tool (SIAN, StructuralIdentifiability.jl, STRIKE-GOLDD) should
  confirm it gold-standard and check for additional symmetry, since this is exactly the less-settled case.
- No surveyed paper fits a transcriptional-vs-translational NAR model to a real dataset with BOTH mRNA and
  protein co-measured over the same time course. Candidate systems: JAK-STAT CIS/SOCS (Bachmann 2011 has
  both), Hes/Her. This combination is a genuine gap and a concrete target (absence in this search is not
  proof of absence; a targeted search is warranted before claiming the gap is fully open).

## Concrete routes forward (ranked)

1. Do the mechanism-discrimination question rigorously and frame it as the result: are M2 and M3
   statistically distinguishable, structurally and practically, from protein-only vs mRNA+protein - using
   the deterministic structural distinguishability (mutual-fit residual vs noise), not the unstable AIC
   rate. This answers the field's open question and is achievable from what we have.
2. Confirm the M2/M3 structural identifiability with a gold-standard tool (StructuralIdentifiability.jl /
   SIAN), upgrading our Fisher+hand-proof and engaging the less-settled autoregulated case.
3. Stretch: connect to a real co-measured mRNA+protein dataset (Bachmann CIS/SOCS) to make the
   discrimination result non-tautological on real data.

## Honest ceiling

Even the best case here is a modest contribution: a clean, correct answer to a narrow open question
(mechanism distinguishability + measurement requirement), not a high-impact result. Sources: Pieschner
2022; Frohlich 2018; Burton 2021; Raue 2009; Villaverde 2019; Massonis-Villaverde-Banga 2023; Gutenkunst
2007; Ryzowicz & Yildirim 2023.
