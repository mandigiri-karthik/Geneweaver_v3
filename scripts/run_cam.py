"""Run the CAM agent on a genome and (optionally) evaluate it against RefSeq.

    python scripts/run_cam.py --genome data/raw/lsakei_23k.fna --ref data/raw/lsakei_23k.gff3

Outputs (in results/<name>/): proteins.faa, genes.gff3, qc_table.tsv,
cam_eval.json, run_log.txt, and a SQLite checkpoint of every graph step.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from geneweaver.cam import compile_cam  # noqa: E402
from geneweaver.evaluation.compare_refseq import compare, format_table  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True, help="genome FASTA")
    ap.add_argument("--ref", default=None, help="reference GFF3 for evaluation (optional)")
    ap.add_argument("--name", default=None, help="run name (default: genome file stem)")
    ap.add_argument("--locus-prefix", default="GW")
    ap.add_argument("--min-protein-len", type=int, default=60)
    a = ap.parse_args(argv)

    name = a.name or Path(a.genome).stem
    outdir = ROOT / "results" / name
    app = compile_cam(checkpoint_db=ROOT / "results" / "checkpoints.sqlite")
    thread_id = f"{name}-{uuid.uuid4().hex[:8]}"
    config = {"configurable": {"thread_id": thread_id}}

    t0 = time.perf_counter()
    final = app.invoke(
        {
            "genome_path": str(Path(a.genome).resolve()),
            "output_dir": str(outdir),
            "locus_prefix": a.locus_prefix,
            "min_protein_len": a.min_protein_len,
            "log": [],
        },
        config=config,
    )
    elapsed = time.perf_counter() - t0

    print("=== CAM run log ===")
    for line in final["log"]:
        print(" ", line)
    print(f"  total: {elapsed:.1f}s   checkpoint thread: {thread_id}")
    print("\n=== CAM stats ===")
    for k, v in final["stats"].items():
        print(f"  {k}: {v}")

    (outdir / "run_log.txt").write_text(
        "\n".join(final["log"]) + f"\nelapsed_s={elapsed:.2f}\nthread_id={thread_id}\n"
    )

    if a.ref:
        metrics = compare(final["outputs"]["gff3"], a.ref)
        metrics["runtime_s"] = round(elapsed, 2)
        (outdir / "cam_eval.json").write_text(json.dumps(metrics, indent=2))
        print("\n=== Evaluation vs reference ===")
        print(format_table(metrics))
        print(f"\nFull metrics: {outdir / 'cam_eval.json'}")


if __name__ == "__main__":
    main()
