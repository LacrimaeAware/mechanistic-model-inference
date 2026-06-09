# Literature and context

## Is this work experimental or mathematical modeling?

The relevant body of work is mathematical modeling and simulation: ordinary and delay differential
equation models of gene and protein dynamics, analyzed for stability, bifurcation, oscillation, and
noise, and compared against experimental data drawn from the literature rather than generated in a wet
lab. The object of this project is therefore models and their statistical properties, not new
experiments.

## First target (recent work, prioritized)

Ryzowicz and Yildirim (2023), "Differential roles of transcriptional and translational negative
autoregulations in protein dynamics," Molecular Omics 19(1):60-71. DOI 10.1039/D2MO00222A; PubMed
36399028.

- Three ODE models of mRNA and protein under an input signal that up-regulates a protein, which then
  represses its own production: (i) no regulation, (ii) transcriptional negative autoregulation,
  (iii) translational negative autoregulation.
- Reported findings: negative autoregulation gives faster responses and quicker return to baseline;
  the transcriptional model is the only one producing oscillations; stochastic simulation orders the
  noise as transcriptional (noisiest), then simplistic, then translational (least).
- Why first: the no-regulation case is the textbook mRNA-to-protein cascade the carried-over tooling
  already validates on, and the three nested models invite an identifiability and model-discrimination
  analysis.
- Status: the author manuscript (31-page accepted version) is held locally and not committed
  (data/raw, gitignored). A first-pass transcription of the three models, the Hill forms, and the
  parameter table is complete. A few equation details still await visual confirmation against the PDF
  before Phase 1 locks: how the input amplitude enters the two regulated models, and the
  nondimensionalization exponents.

## Recent work by the first target's author (candidate later targets)

- Yildirim (2025), "Cellular Decision Making and Signal Adaptation Via Incoherent Feedforward and
  Negative Feedback Loops," Differential Equations and Dynamical Systems. DOI
  10.1007/s12591-025-00735-z. Two three-state ODE motifs (an incoherent feedforward loop and a
  negative feedback loop) that both produce adaptation but differ in signature: linear versus
  saturating dose response, exponential versus oscillatory return. Held locally (data/raw, gitignored).
  Cites the first target. Strong candidate for the Phase 5 second target: it concretizes the
  cellular-adaptation placeholder in the roadmap and mirrors the first target's two-mechanism
  discrimination question.
- Yildirim, Brew and Ay (2025), "Regulatory Effects of Cooperativity and Signal Profile on Adaptive
  Dynamics in Incoherent Feedforward Loop Networks," In Silico Biology 16(1). DOI
  10.1177/14343207241306092. Companion on the effect of cooperativity (Hill coefficient) and signal
  profile on feedforward adaptation. Not held locally. Identifiers verified via Crossref.

## Additional source (textbook chapter)

Robeva and Yildirim (2024), Chapter 2, "Bistability in bacterial genetic networks: a comparison of
differential equations and Boolean network models," in Mathematical Concepts and Methods in Modern
Biology (2nd edition, ISBN 9780443296529), pp. 47-95. DOI 10.1016/B978-0-443-29652-9.00001-4. A
comparison of continuous ODE and discrete Boolean models of bistability in the lactose and arabinose
operons. Relevant as context for model choice and as a source of candidate targets. The minimal
three-variable arabinose ODE model (delay-free, bistable in external arabinose) is the most tractable
identifiability target it contains; the lactose models are delay-differential, which complicates
standard identifiability tooling. Full parameter tables are not reproduced in the chapter (they are in
the cited operon sources). Status: local copy extracted and summarized; bistability windows and
steady-state tables recorded for later use.

## Related modeling work on autoregulation (independent of the first target's author)

- Murray (2023), "Autoregulation of Transcription and Translation: A Qualitative Analysis," Bulletin
  of Mathematical Biology 85(7):57. DOI 10.1007/s11538-023-01143-6; PubMed 37233955.
- Murray, Ocana, Meijer and Dale (2021), "Auto-Regulation of Transcription and Translation:
  Oscillations, Excitability and Intermittency," Biomolecules 11(11):1566. DOI 10.3390/biom11111566;
  PMC8615617.

Both are independent works, not by the first target's author, on the same model family: a two-variable
mRNA/protein ODE system with Hill-type transcriptional self-repression and nonlinear regulation of
translation, producing oscillation, excitability, and homeostasis. The 2021 paper models the Hes/Her
system and is the precursor to the 2023 analysis. They share the first target's modeling form and are
natural comparators for the model-discrimination question. Identifiers verified via Crossref and
PubMed.

## Carried-over context (verified in the sibling project; 23 of 25 claims confirmed there)

Mechanistic-modeling corpus:
- Yildirim and Mackey (2003), Biophysical Journal 84:2841-2851: a delay-differential-equation model of
  the lac operon with bistability via a cusp bifurcation.
- Dyjack, Azeredo-Tseng and Yildirim (2017), Molecular BioSystems 13(7):1323-1335: cellular adaptation
  in the yeast MAPK/pheromone response via two negative-feedback phosphatases.
- Violin et al. (2008), JBC 283(5):2949-2961: beta2-adrenergic receptor desensitization (adaptation).

Cellular-adaptation theory (control-theoretic backbone):
- Sontag (2003/2004), Systems and Control Letters: the internal model principle.
- Robust perfect adaptation as integral feedback; antithetic integral feedback (Briat, Gupta and
  Khammash 2016, Cell Systems; Aoki et al. 2019, Nature 570:533-537).
- Araujo and Liotta (2018), Nature Communications 9:1757: topological classification of robust
  perfect adaptation networks into opposer and balancer modules.
- Gupta and Khammash (2022), PNAS 119:e2207802119: robust-perfect-adaptation constraints differ
  between deterministic and stochastic descriptions.

The open gap (the project's reason to exist): no rigorous parameter-identifiability and inference
analysis of these small adaptation and gene-expression models has been published. The dynamics have
been built and analyzed; the estimability questions (which parameters can be recovered, how well, from
what measurements) have not. That gap is the statistics and quantitative contribution.

## Privacy note

The citations above are ordinary scholarly reference. This project does not represent collaboration
with, mentorship by, or endorsement from any author, and does not state any author's employer or role.
See AGENTS.md.
