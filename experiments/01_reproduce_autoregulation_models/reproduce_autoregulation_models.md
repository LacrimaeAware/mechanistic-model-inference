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

Not started. The no-regulation baseline (mRNA to protein) already exists as `simulate_mrna_protein`.
The two regulation models are stubbed in `models.py`, pending transcription from the paper.

## Notes

- Do not reconstruct the equations from memory; use the paper text.
- The purpose of reproduction is a trusted forward model, which is the prerequisite for the
  identifiability and inference work in later phases.
