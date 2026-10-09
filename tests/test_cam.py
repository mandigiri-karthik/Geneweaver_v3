"""Tests for the CAM agent and the RefSeq evaluation.

A synthetic ~60 kb genome with planted genes on both strands gives a known
ground truth without needing internet access.
"""
from __future__ import annotations

import random
import sys
from pathlib import Path

import pytest
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from geneweaver.cam import compile_cam  # noqa: E402
from geneweaver.cam import nodes  # noqa: E402
from geneweaver.evaluation.compare_refseq import compare  # noqa: E402

STOPS = {"TAA", "TAG", "TGA"}
SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS]


def _orf(rng: random.Random, n_codons: int) -> str:
    return "ATG" + "".join(rng.choice(SENSE) for _ in range(n_codons)) + "TAA"


def _intergenic(rng: random.Random, n: int) -> str:
    # AT-rich spacer so it is unlikely to contain long ORFs
    return "".join(rng.choices("ACGT", weights=[35, 15, 15, 35], k=n))


@pytest.fixture(scope="module")
def synthetic(tmp_path_factory):
    rng = random.Random(42)
    d = tmp_path_factory.mktemp("synthetic")
    seq, genes, pos = [], [], 0

    def add(s):
        nonlocal pos
        seq.append(s)
        pos += len(s)

    add(_intergenic(rng, 300))
    for i in range(60):
        add(_intergenic(rng, rng.randint(80, 200)) + "AGGAGG" + _intergenic(rng, 8))
        orf = _orf(rng, rng.randint(150, 400))
        strand = "+" if i % 3 else "-"
        start = pos + 1
        add(orf if strand == "+" else str(Seq(orf).reverse_complement()))
        genes.append((start, pos, strand))
    add(_intergenic(rng, 300))
    genome = "".join(seq)

    fna = d / "synth.fna"
    fna.write_text(">synth\n" + "\n".join(genome[i:i + 70] for i in range(0, len(genome), 70)) + "\n")
    ref = d / "synth_ref.gff3"
    with ref.open("w") as fh:
        fh.write("##gff-version 3\n")
        for n, (s, e, st) in enumerate(genes, 1):
            fh.write(f"synth\tplanted\tCDS\t{s}\t{e}\t.\t{st}\t0\tID=cds{n};locus_tag=REF_{n}\n")
    return {"fna": fna, "ref": ref, "n": len(genes), "dir": d}


def test_graph_runs_end_to_end(synthetic):
    app = compile_cam()
    out = synthetic["dir"] / "out"
    final = app.invoke({"genome_path": str(synthetic["fna"]), "output_dir": str(out), "log": []})
    assert final["stats"]["n_genes"] > 0
    for key in ("proteins", "gff3", "qc_table"):
        assert Path(final["outputs"][key]).exists()
    # every node wrote to the audit log
    assert [l.split(":")[0] for l in final["log"]][:4] == ["load_genome", "call_genes", "translate", "qc_check"]
    # Biopython and Pyrodigal agree on complete genes
    assert not any("translation_mismatch" in g["qc_flags"] for g in final["genes"])


def test_recovers_planted_genes(synthetic):
    app = compile_cam()
    out = synthetic["dir"] / "eval"
    final = app.invoke({"genome_path": str(synthetic["fna"]), "output_dir": str(out), "log": []})
    m = compare(final["outputs"]["gff3"], synthetic["ref"])
    assert m["n_reference_cds"] == synthetic["n"]
    assert m["sensitivity_pct"] >= 90


def test_checkpointing(synthetic, tmp_path):
    app = compile_cam(checkpoint_db=tmp_path / "ckpt.sqlite")
    cfg = {"configurable": {"thread_id": "t1"}}
    app.invoke({"genome_path": str(synthetic["fna"]), "output_dir": str(tmp_path / "o"), "log": []}, cfg)
    history = list(app.get_state_history(cfg))
    assert len(history) >= 6  # one snapshot per step


def test_qc_flags_short_and_internal_stop():
    state = {"min_protein_len": 60, "genes": [
        {"locus_tag": "A", "protein": "M" * 100, "partial_begin": False, "partial_end": False,
         "start_codon": "ATG", "rbs_motif": "GGAG", "qc_flags": []},
        {"locus_tag": "B", "protein": "M" * 20, "partial_begin": False, "partial_end": False,
         "start_codon": "GTG", "rbs_motif": None, "qc_flags": []},
        {"locus_tag": "C", "protein": "MK*KL" * 30, "partial_begin": False, "partial_end": True,
         "start_codon": "ATG", "rbs_motif": "GGAG", "qc_flags": []},
    ]}
    upd = nodes.qc_check(state)
    flags = {g["locus_tag"]: set(g["qc_flags"]) for g in upd["genes"]}
    assert flags["A"] == set()
    assert {"short_protein", "non_atg_start", "no_rbs_motif"} <= flags["B"]
    assert {"internal_stop", "partial"} <= flags["C"]
    assert upd["flagged_for_review"] == ["B", "C"]


def test_compare_exact_vs_start_error(tmp_path):
    ref = tmp_path / "ref.gff3"
    pred = tmp_path / "pred.gff3"
    ref.write_text(
        "c\tr\tCDS\t100\t400\t.\t+\t0\tID=r1;locus_tag=R1\n"
        "c\tr\tCDS\t1000\t1300\t.\t-\t0\tID=r2;locus_tag=R2\n"
        "c\tr\tCDS\t2000\t2300\t.\t+\t0\tID=r3;locus_tag=R3\n"
        "c\tr\tCDS\t3000\t3300\t.\t+\t0\tID=r4;pseudo=true\n"
    )
    pred.write_text(
        "c\tp\tCDS\t100\t400\t.\t+\t0\tID=p1\n"     # exact
        "c\tp\tCDS\t1000\t1330\t.\t-\t0\tID=p2\n"   # same stop (1000), start 30 bp upstream
        "c\tp\tCDS\t5000\t5300\t.\t+\t0\tID=p3\n"   # extra
    )
    m = compare(pred, ref)
    assert (m["exact_match"], m["start_errors"], m["missed"], m["extra"]) == (1, 1, 1, 1)
    assert m["n_reference_pseudo_skipped"] == 1
    assert m["start_shift_longer_vs_shorter"]["cam_longer"] == 1
