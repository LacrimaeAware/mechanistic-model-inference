# Experiment 01: reproduce the negative-autoregulation models

## Goal

Transcribe and reproduce the three ODE models from the 2023 negative-autoregulation paper
(docs/literature.md): (i) no regulation, (ii) transcriptional negative autoregulation, (iii)
translational negative autoregulation. Match the published deterministic behavior before any
identifiability or inference claim is made.

## Plan

1. Transcribe the exact right-hand sides, Hill exponents, and parameter values from the paper text
   (not from memory) into `mechanistic_inference.models`.
2. Simulate each model under the paper's input signal and reproduce the reported behavior: faster
   response and quicker return to baseline under regulation, and oscillation in the transcriptional
   model only.
3. Optional: a Gillespie stochastic simulation to reproduce the noise ordering (transcriptional
   noisiest, then simplistic, then translational).
4. Save figures to `results/` and record the numbers in this note.

## Status

Deterministic reproduction done. The three dimensionless models (Eqs. 4-7) are implemented in
`src/mechanistic_inference/models.py`, transcribed from the paper and confirmed against the rendered
PDF pages (the equations are not committed; the manuscript stays in `data/raw`, gitignored). Two
ambiguities in the extracted text were resolved by reading the PDF: the Hill term is `1/(1+(kappa P)^n)`
(Eq. 3), and the input signal scales transcription multiplicatively as `(1+S_a) F(p)` (S_a in [1,4] is
a 2- to 5-fold increase, p. 15-16). The optional Gillespie stochastic reproduction is also done.

Regenerate the deterministic results with `python experiments/01_reproduce_autoregulation_models/reproduce.py`
(writes `results/fig01_oscillation.png`, `results/fig02_response_speed.png`,
`results/reproduce_01_summary.txt`) and the stochastic results with `reproduce_stochastic.py` (writes
`results/fig03_stochastic_cv.png`, `results/reproduce_01_stochastic_summary.txt`). Correctness is
covered by `tests/test_models.py` and `tests/test_stochastic.py` (full suite 17 passing).

## Results (regenerated)

- M1 numerical integration matches the closed form (Eq. 8): max abs error 3.0e-8.
- Oscillation only in the transcriptional model (theta=0.6, n=8): M2 discriminant -12.74 at
  1/kappa_tilde_M=0.5 (damped oscillation, protein overshoot 0.107) and +0.16 at 1/kappa_tilde_M=6
  (monotone). M1 and M3 overshoot ~1e-13 (never oscillate). Reproduces Fig. 4.
- Faster response and quicker return under regulation (step S_a=3, theta=0.0571, 1/kappa=50 molecules):
  half-response time tau_h = 13.2 (M1), 3.6 (M2), 4.1 (M3); return time tau_r = 48.4 (M1), 9.2 (M2),
  11.4 (M3). M1 is slowest by ~3.6x (response) and ~5x (return); M3 is slightly slower than M2.
  Reproduces the Fig. 5 ordering.
- Two structural facts recorded for later identifiability work: M3 mRNA is identical to M1 mRNA over
  the step (max diff 6.2e-9), so mRNA cannot separate M1 from M3; and M2 and M3 share the same protein
  steady state under symmetric feedback (kappa_tilde_M=kappa_tilde_P, n_M=n_P), a discrimination
  degeneracy in the protein-level observable.
- Stochastic noise ordering reproduced (exact SSA, 300 runs, 1/kappa=300 molecules, S_a=0): protein
  CV M2=0.153 > M1=0.141 > M3=0.091 (paper: 0.16, 0.14, 0.10); mRNA CV M2=0.96 > M1=0.60 ~ M3=0.58
  (paper: 0.89, 0.56, 0.56). The ordering and the protein CVs match within ensemble noise. Caveat:
  the Hill-in-propensity formulation is a modeling choice the paper does not write out, so the exact
  CVs are not claimed identical to theirs; M2 mRNA CV runs somewhat higher than reported.

## Open items before Phase 2

- A few equation details were confirmed visually but the Jacobian/eigenvalue algebra (Eq. 9, h1/h2,
  pp. 10-12) was used only as an oscillation check, not re-derived; revisit if needed.
- The stochastic propensity formulation for the Hill repression is one of several plausible choices;
  an alternative (promoter-state or occupancy-based) could shift the CV values.

## Notes

- Do not reconstruct the equations from memory; use the paper text.
- The purpose of reproduction is a trusted forward model, which is the prerequisite for the
  identifiability and inference work in later phases.
