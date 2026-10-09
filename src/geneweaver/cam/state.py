"""Typed state for the CAM (Contig Analysis Manager) agent.

CAM is the heir of AGMIAL's Contig Analysis Manager (Bryson et al., NAR 2006):
it takes a genome and produces gene models and translated proteins, which are
handed to the PAM agent for functional annotation.
"""
from __future__ import annotations

import operator
from typing import Annotated, TypedDict


class GeneModel(TypedDict):
    locus_tag: str
    contig: str
    start: int          # 1-based, inclusive (GFF convention)
    end: int            # 1-based, inclusive
    strand: str         # "+" or "-"
    start_codon: str    # e.g. "ATG", "GTG", "TTG"
    partial_begin: bool
    partial_end: bool
    rbs_motif: str | None
    score: float
    protein: str        # translated sequence, without the trailing stop "*"
    qc_flags: list[str]  # empty list = passed QC


class CAMState(TypedDict, total=False):
    # inputs
    genome_path: str
    output_dir: str
    locus_prefix: str
    min_protein_len: int

    # filled in by nodes
    contigs: dict[str, str]          # contig id -> DNA sequence
    genes: list[GeneModel]
    flagged_for_review: list[str]    # locus tags sent to PAM / human review
    outputs: dict[str, str]          # name -> file path
    stats: dict[str, float]

    # append-only audit trail: every node writes one line (evidence of process)
    log: Annotated[list[str], operator.add]
