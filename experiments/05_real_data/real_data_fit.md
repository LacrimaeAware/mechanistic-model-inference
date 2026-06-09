# Experiment 05: fit a small mRNA->protein model to real single-cell data

## Goal

Fit a small mRNA->protein ODE to real single-cell GFP translation data and report which parameters are
identifiable, then compare the answer to a published structural-identifiability analysis of the same
data. This is the step that takes the project off synthetic data: a real fit, against a known benchmark,
on exactly the protein-observed / mRNA-latent identifiability question the project studies.

## Data

Frohlich et al. (2018), "Multi-experiment nonlinear mixed effect modeling of single-cell translation
kinetics after transfection," npj Systems Biology and Applications, DOI 10.1038/s41540-018-0079-7.
Chemically modified mRNA encoding eGFP is transfected as a bolus at t=0; GFP fluorescence is imaged in
single cells while the mRNA (latent) degrades. Source: Zenodo 10.5281/zenodo.1228899 (code.zip),
CC-BY-SA-4.0; the raw zip is in data/raw (gitignored), downloaded on demand by
`experiments/05_real_data/load_and_inspect.py`.

Loaded and inspected (file 20160427_mean_eGFP.xlsx): 180 timepoints every 600 s over 30 h, 1079
single-cell traces. The mean trace rises from ~8 a.u. to a peak of ~2489 a.u. near 14.7 h, then decays,
which is the expected translation signature (the mRNA bolus depletes, so protein production wanes and
protein decays). Figure: `results/fig07_frohlich_data.png`. Cell-to-cell variability is large (single
cells peak from a few hundred to ~6500 a.u.), so single-cell fits and a mean-trace fit will both be run.

## Model to implement

The translation model used for this data (Leonhardt et al. 2014 / Frohlich et al. 2018): an mRNA bolus
with a translation onset time, translation into protein, and protein decay, observed through a
fluorescence scale and offset. Schematically, for t >= t0:

    dm/dt = -delta * m ,   m(t0) = m0
    dp/dt = k * m - beta * p ,   p(t0) = 0
    observed GFP(t) = scale * p(t) + offset

Known identifiability structure (the reason this dataset is a good benchmark): m0, the translation rate
k, and the fluorescence scale enter the observed signal only through their product, so they are not
separately identifiable from protein (GFP) alone, only the combination is. This is the same
transcription/translation-rate degeneracy theme as the rest of the project, in a real model with a
published answer (Pieschner, Hasenauer & Fuchs 2022, "Identifiability analysis for models of the
translation kinetics after mRNA transfection," PMC9110294).

## Plan

1. Implement the translation model in `mechanistic_inference.models` (with the onset time and the
   scale/offset observation), as a `simulate(theta, t)` closure for the identifiability tooling.
2. Fit it to the mean trace and to a sample of single-cell traces with `fit_least_squares`; report the
   fit quality and the recovered parameters.
3. Run Fisher rank and profile likelihood on the real fit; name the non-identifiable combination and
   check it against the published m0 * k * scale degeneracy.
4. Deliverable: a figure of the fit over the data and an identifiability statement for the real model,
   stated as a comparison to the published benchmark.

## Status

Data downloaded, loaded, and visualized (this note; `load_and_inspect.py`; `fig07_frohlich_data.png`).
The model and the fit are not yet implemented. This is the immediate next step.

## Notes

- The benchmark gives a correctness check the synthetic work could not: if the pipeline does not recover
  the published degeneracy on this real data, the pipeline (or its application) is wrong.
- The Frohlich model differs from the project's M1 in that mRNA is a decaying bolus (no production term)
  and the observation has a scale and offset; the identifiability question (protein-only) is the same.
