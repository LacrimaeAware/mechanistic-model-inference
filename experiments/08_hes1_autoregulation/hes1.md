# Experiment 08: real Hes1 autoregulation data - characterization and an honest model-fit assessment

## What this is

Hes1 is a real negative-autoregulation gene (it represses its own transcription), in the family of the
source paper's transcriptional model. The goal was to do the "fit real data, recover parameters"
deliverable on an actual autoregulation gene. Data: GPosc repo (Phillips et al. 2017),
`data/raw/hes1/Hes1_example.csv` (gitignored).

## What the data actually is (looked at, not assumed)

12 single-cell fluorescent-reporter traces, 0.5 h sampling, each tracked for a different duration (median
~37 h, range 24-46 h; the column goes NaN when the cell is lost, so each cell is handled on its own valid
window). Background-subtracted, the traces start high (~300-600 a.u.) and decay over the first hours to
~50-150 a.u., then fluctuate irregularly. A crude autocorrelation gives an ~11 h dominant timescale, but
that is NOT trustworthy: the traces are noisy and the source paper detected Hes1 oscillations only with
Gaussian-process methods, precisely because simple period estimates fail here. Figure
`results/fig12_hes1_data.png`.

## Honest assessment: the current models do not fit this data, and forcing a fit would mislead

- Shape: the Hes1 traces DECAY from an initial high; the project's mRNA->protein autoregulation models
  RISE from zero to a steady state. Opposite qualitative shape.
- Oscillation: sustained Hes1 oscillation requires a transcriptional DELAY (Monk 2003); the project's
  models are non-delay and can only damp, so they cannot produce the ongoing fluctuations.
- Noise and confounds: the single-cell traces are noisy and likely carry experimental confounds (a
  photobleaching-like initial decay, real biological noise) that the clean models do not represent; the
  source handled this with a stochastic Gaussian-process model.

So a meaningful Hes1 analysis is a substantial build - a delayed autoregulation model plus a
stochastic/experimental-noise treatment - not a quick fit with the current tooling. Fitting M2 or M3 to
these traces would produce numbers without meaning, which is the failure mode this project is trying to
avoid. Recorded as an honest stop, not a result.

## Options from here

- Build a delayed autoregulation model (Monk 2003) and a noise treatment as a deliberate, scoped effort:
  the genuinely relevant direction, but with real risk of a long, possibly inconclusive effort given the
  data's noise (the reason the source used Gaussian processes).
- Treat the Frohlich single-cell parameter recovery (experiment 05) as the project's real-data
  deliverable and consolidate, with Hes1 noted as future work that needs a delay model.
