"""Nodes of the CAM LangGraph subgraph.

Each node is a plain function: it takes the current CAMState and returns a
partial update. Keeping them pure makes them easy to unit-test (tests/test_cam.py).
"""
from __future__ import annotations

from pathlib import Path

import pyrodigal
from Bio import SeqIO
from Bio.Seq import Seq

from .state import CAMState, GeneModel

# Pyrodigal recommends >= 100 kb to self-train; below this we use metagenomic mode.
MIN_TRAINING_LENGTH = 100_000

# QC flags that send a gene to review (PAM evidence check, later a human).
REVIEW_FLAGS = {"internal_stop", "partial", "short_protein", "translation_mismatch"}


def load_genome(state: CAMState) -> dict:
    """Read every contig from a FASTA file."""
    path = Path(state["genome_path"])
    contigs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(path, "fasta")}
    if not contigs:
        raise ValueError(f"No sequences found in {path}")
    total = sum(len(s) for s in contigs.values())
    return {
        "contigs": contigs,
        "log": [f"load_genome: {len(contigs)} contig(s), {total:,} bp from {path.name}"],
    }


def call_genes(state: CAMState) -> dict:
    """Predict protein-coding genes with Pyrodigal (Prodigal, translation table 11)."""
    contigs = state["contigs"]
    total = sum(len(s) for s in contigs.values())
    meta = total < MIN_TRAINING_LENGTH

    finder = pyrodigal.GeneFinder(meta=meta)
    if not meta:
        # Train once on all contigs so the model reflects the whole genome.
        finder.train(*(s.encode() for s in contigs.values()), translation_table=11)

    prefix = state.get("locus_prefix", "GW")
    genes: list[GeneModel] = []
    n = 0
    for contig_id, seq in contigs.items():
        for g in finder.find_genes(seq.encode()):
            n += 1
            genes.append(
                GeneModel(
                    locus_tag=f"{prefix}_{n * 5:05d}",  # NCBI-style spacing of 5
                    contig=contig_id,
                    start=g.begin,
                    end=g.end,
                    strand="+" if g.strand == 1 else "-",
                    start_codon=g.start_type,
                    partial_begin=g.partial_begin,
                    partial_end=g.partial_end,
                    rbs_motif=g.rbs_motif,
                    score=round(g.score, 2),
                    protein=g.translate(translation_table=11).rstrip("*"),
                    qc_flags=[],
                )
            )
    mode = "metagenomic (genome too short to train)" if meta else "single-genome (self-trained)"
    return {"genes": genes, "log": [f"call_genes: {len(genes)} CDS predicted, mode={mode}"]}


def translate(state: CAMState) -> dict:
    """Independently re-translate each CDS with Biopython (table 11).

    Pyrodigal already returns a protein; translating again with a second
    library is a cheap consistency check. Disagreements are flagged, not fixed.
    """
    contigs = state["contigs"]
    genes = [dict(g) for g in state["genes"]]
    mismatches = 0
    for g in genes:
        if g["partial_begin"] or g["partial_end"]:
            continue  # frame at a contig edge is not reliably translatable
        dna = Seq(contigs[g["contig"]][g["start"] - 1 : g["end"]])
        if g["strand"] == "-":
            dna = dna.reverse_complement()
        try:
            prot = str(dna.translate(table=11, cds=True))
        except Exception:  # Bio.Data.CodonTable.TranslationError and friends
            prot = None
        if prot != g["protein"]:
            g["qc_flags"] = g["qc_flags"] + ["translation_mismatch"]
            mismatches += 1
    return {
        "genes": genes,
        "log": [f"translate: Biopython cross-check, {mismatches} mismatch(es)"],
    }


def qc_check(state: CAMState) -> dict:
    """Rule-based QC. Flags genes whose model looks suspect for PAM review."""
    min_len = state.get("min_protein_len", 60)
    genes = [dict(g) for g in state["genes"]]
    for g in genes:
        flags = list(g["qc_flags"])
        if "*" in g["protein"]:
            flags.append("internal_stop")
        if g["partial_begin"] or g["partial_end"]:
            flags.append("partial")
        if len(g["protein"]) < min_len:
            flags.append("short_protein")
        if g["start_codon"] != "ATG":
            flags.append("non_atg_start")  # informational: GTG/TTG are legitimate
        if not g["rbs_motif"] or g["rbs_motif"] == "None":
            flags.append("no_rbs_motif")    # informational
        g["qc_flags"] = sorted(set(flags))

    flagged = [g["locus_tag"] for g in genes if REVIEW_FLAGS & set(g["qc_flags"])]
    return {
        "genes": genes,
        "flagged_for_review": flagged,
        "log": [f"qc_check: {len(flagged)} of {len(genes)} gene(s) flagged for review"],
    }


def route_after_qc(state: CAMState) -> str:
    """Conditional edge. Today both branches end at export; later the
    'review' branch hands flagged genes to the PAM agent (the CAM <-> PAM loop)."""
    return "review" if state.get("flagged_for_review") else "clean"


def mark_for_review(state: CAMState) -> dict:
    """Placeholder for the CAM -> PAM hand-off of suspect gene models."""
    return {"log": [f"mark_for_review: {len(state['flagged_for_review'])} gene(s) queued for PAM"]}


def export(state: CAMState) -> dict:
    """Write proteins (FASTA), gene models (GFF3) and a per-gene QC table (TSV)."""
    out = Path(state["output_dir"])
    out.mkdir(parents=True, exist_ok=True)
    genes = state["genes"]

    faa = out / "proteins.faa"
    with faa.open("w") as fh:
        for g in genes:
            fh.write(f">{g['locus_tag']} {g['contig']}:{g['start']}-{g['end']}({g['strand']})\n")
            for i in range(0, len(g["protein"]), 60):
                fh.write(g["protein"][i : i + 60] + "\n")

    gff = out / "genes.gff3"
    with gff.open("w") as fh:
        fh.write("##gff-version 3\n")
        for cid, seq in state["contigs"].items():
            fh.write(f"##sequence-region {cid} 1 {len(seq)}\n")
        for g in genes:
            partial = ""
            if g["partial_begin"] or g["partial_end"]:
                partial = f";partial={int(g['partial_begin'])}{int(g['partial_end'])}"
            attrs = (
                f"ID={g['locus_tag']};locus_tag={g['locus_tag']};start_codon={g['start_codon']}"
                f";rbs_motif={g['rbs_motif']}{partial}"
            )
            fh.write(
                f"{g['contig']}\tGeneWeaver-CAM\tCDS\t{g['start']}\t{g['end']}\t"
                f"{g['score']}\t{g['strand']}\t0\t{attrs}\n"
            )

    tsv = out / "qc_table.tsv"
    with tsv.open("w") as fh:
        fh.write("locus_tag\tcontig\tstart\tend\tstrand\tlength_aa\tstart_codon\tscore\tqc_flags\n")
        for g in genes:
            fh.write(
                f"{g['locus_tag']}\t{g['contig']}\t{g['start']}\t{g['end']}\t{g['strand']}\t"
                f"{len(g['protein'])}\t{g['start_codon']}\t{g['score']}\t{','.join(g['qc_flags'])}\n"
            )

    lengths = [len(g["protein"]) for g in genes] or [0]
    stats = {
        "n_genes": len(genes),
        "n_flagged": len(state.get("flagged_for_review", [])),
        "mean_protein_len": round(sum(lengths) / len(lengths), 1),
        "pct_atg_start": round(100 * sum(g["start_codon"] == "ATG" for g in genes) / max(len(genes), 1), 1),
    }
    return {
        "outputs": {"proteins": str(faa), "gff3": str(gff), "qc_table": str(tsv)},
        "stats": stats,
        "log": [f"export: wrote {faa.name}, {gff.name}, {tsv.name} to {out}"],
    }
