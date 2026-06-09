# Parameter identifiability and inference for mechanistic models of gene and protein dynamics

A study of a specific statistical question about small mechanistic biology models: given a published
ODE model of gene or protein dynamics and realistic data, which parameters can actually be estimated,
how well, and from what measurements. The method is the statistician's toolkit (maximum likelihood,
Fisher information, profile likelihood, and later Bayesian inference), applied to models that are
usually built and analyzed dynamically but rarely analyzed for parameter identifiability.

This is the early-stage start of a project that is expected to grow from one concrete target outward,
the same way its sibling project did. The goal is a modest, correct, useful result, not a high-impact
claim.

## In plain terms

A mechanistic model is a set of equations for how molecules change over time (for example, how much
messenger RNA and protein a gene produces, and how the protein feeds back to slow its own production).
Researchers build these models and study their behavior (oscillations, stability). A separate, less
asked question is statistical: if you only have noisy measurements, can you actually pin down the
numbers in the equations? Often you cannot, and the interesting result is showing exactly which
numbers are recoverable, which are not, and what extra measurement would fix that. That is the
question this project works on.

## First target

The first concrete model is a recent negative-autoregulation paper that builds three nested ODE models
of an mRNA and protein under an input signal: (i) no regulation, (ii) transcriptional negative
autoregulation, and (iii) translational negative autoregulation, and compares their response speed,
oscillation, and noise. Full citation and summary in [`docs/literature.md`](docs/literature.md). It is
a good first target because the no-regulation case is the textbook mRNA-to-protein cascade the carried
over tooling is already validated on, and the three models invite an identifiability and
model-discrimination analysis: can the two regulation modes be told apart from data, and under what
measurement scheme.

## Goals

Learning goals (the reason for starting here): ODE simulation, what each parameter means physically,
least-squares and maximum-likelihood fitting, Fisher information, profile likelihood, and later
Bayesian / MCMC inference. Project goals: reproduce a published model from its equations, report which
parameters are identifiable from which measurements, and attempt a model-discrimination result. Detail
and the realistic-versus-stretch split are in [`docs/goals.md`](docs/goals.md) and
[`docs/roadmap.md`](docs/roadmap.md).

## Status

Scaffold. The identifiability pipeline (simulate, MLE, profile likelihood, Fisher information) is
carried over and validated on the textbook mRNA-to-protein model (`tests/test_identifiability.py`,
213-test sibling suite reduced here to the relevant tooling). No target model is reproduced yet; that
is the first experiment.

## Layout

```text
mechanistic-model-inference/
├── src/mechanistic_inference/   # identifiability/inference tooling (validated) and model definitions
├── experiments/                 # one folder per experiment, each with a write-up and a script
├── docs/                        # goals, research questions, roadmap, literature, conventions
└── tests/                       # correctness tests for the tooling
```

## Relation to prior work

Sibling to `stable-grn-inference`, which mapped the limits of directed gene-network edge recovery and
concluded that direction is not recoverable from static data and simple baselines are not beaten. This
project takes the forward direction identified there: not data-driven network inference, but
parameter identifiability and inference for individual mechanistic models, where the bar is a correct
and useful analysis rather than beating a benchmark.

## Citing and privacy

Public-facing. Published papers are cited by author and year in [`docs/literature.md`](docs/literature.md),
which is normal scholarly reference. The project does not represent any collaboration, affiliation, or
endorsement, and does not state any individual's employer or role. See [`AGENTS.md`](AGENTS.md).
