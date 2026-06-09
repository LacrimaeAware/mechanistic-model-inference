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
- Status: abstract and bibliographic details confirmed by web search. The exact equations, Hill
  exponents, and parameter values must be transcribed from the paper text before Phase 1.

## Additional source (textbook chapter)

Chapter 2, "Stability in bacterial genetic networks: a comparison of differential equations and
Boolean network models," in Mathematical Concepts and Methods in Modern Biology (2nd edition, 2024,
ISBN 9780443296529). A direct comparison of continuous ODE and discrete Boolean models for the same
genetic-network stability question. Relevant as context for model choice and as a possible second
target. Status: a local copy exists; it needs text extraction and has not yet been summarized from the
source.

## Related modeling work on autoregulation (found, to verify against the sources)

- "Autoregulation of Transcription and Translation: A Qualitative Analysis," Bulletin of Mathematical
  Biology (2023).
- "Auto-Regulation of Transcription and Translation: Oscillations, Excitability and Intermittency"
  (PMC8615617).

These overlap the first target's theme and may inform the model forms and the oscillation analysis.

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
