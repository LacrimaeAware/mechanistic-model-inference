"""Mechanistic ODE models for the project's targets.

The no-regulation mRNA -> protein baseline already exists as ``simulate_mrna_protein`` in
:mod:`mechanistic_inference.identifiability`. The two negative-autoregulation models below are the
first target (the 2023 negative-autoregulation paper; full citation in docs/literature.md): an input
signal up-regulates a protein, which then represses its own production either transcriptionally or
translationally.

These are left unimplemented on purpose. The exact right-hand sides, Hill exponents, and parameter
values must be transcribed from the paper (not reconstructed from memory) and matched against its
published figures before use. See experiments/01_reproduce_autoregulation_models.
"""
from __future__ import annotations


def simulate_transcriptional_nar(*args, **kwargs):
    """Model (ii): transcriptional negative autoregulation (repression of mRNA production).

    TODO: implement from the published equations and validate against the paper's figures.
    """
    raise NotImplementedError(
        "Transcribe the equations and parameters from the paper first; see docs/literature.md."
    )


def simulate_translational_nar(*args, **kwargs):
    """Model (iii): translational negative autoregulation (repression of protein production from mRNA).

    TODO: implement from the published equations and validate against the paper's figures.
    """
    raise NotImplementedError(
        "Transcribe the equations and parameters from the paper first; see docs/literature.md."
    )
