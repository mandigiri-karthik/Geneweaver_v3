"""Sanity checks for score_gff.py. Exits with status 1 if any check fails.

Usage (from Phase-0/):
    python3 sanity_checks.py [--outdir sanity]

1. Self-score: the reference GFF3 with its UniProt EC numbers added (ec_number=) is scored against
   the reference. Every recovery and agreement metric must be 100%, every error metric 0%.
2. Deletion: 10 random reference genes (seed 42) are removed from a copy of that file. The scorer must
   report exactly those 10 genes as missed, and nothing else must change.

Uses only the Python standard library.
"""
import argparse
import csv
import random
import sys
from pathlib import Path

import score_gff

DELETE_COUNT = 10
DELETE_SEED = 42

MUST_BE_ONE = ["recall", "precision", "start_exact", "symbol_match", "ec_match", "ec_set_equal",
               "product_match", "agrees"]
MUST_BE_ZERO = ["false_positives", "on_reference_pseudogenes", "hypothetical_where_known", "wrong_function"]


def write_reference_with_ec(reference_gff, reference_function, output_path, drop_ids=frozenset()):
    """Copy the reference GFF3, adding ec_number= to each CDS row and dropping CDS rows in drop_ids."""
    with open(reference_function) as fh:
        ecs = {row["cds_id"]: row["ec_numbers"] for row in csv.DictReader(fh, delimiter="\t")}

    with open(reference_gff) as src, open(output_path, "w") as out:
        for line in src:
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(fields) != 9 or fields[2] != "CDS":
                out.write(line)
                continue
            cds_id = score_gff.parse_attributes(fields[8])["ID"]
            if cds_id in drop_ids:
                continue
            if ecs.get(cds_id):
                fields[8] = fields[8].rstrip(";") + ";ec_number=" + ecs[cds_id].replace(";", ",")
            out.write("\t".join(fields) + "\n")


def check(label, passed, detail=""):
    print(f"  [{'PASS' if passed else 'FAIL'}] {label}{(': ' + detail) if detail else ''}")
    return passed


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--reference-gff", type=Path, default=score_gff.DEFAULT_REFERENCE_GFF)
    parser.add_argument("--reference-function", type=Path, default=score_gff.DEFAULT_REFERENCE_FUNCTION)
    parser.add_argument("--outdir", type=Path, default=score_gff.HERE / "sanity")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    genes, pseudo_spans = score_gff.load_reference(args.reference_gff, args.reference_function)
    all_passed = True

    # 1. Self-score
    print("Check 1: reference (with EC) scored against itself")
    self_gff = args.outdir / "reference_with_ec.gff3"
    write_reference_with_ec(args.reference_gff, args.reference_function, self_gff)
    metrics = score_gff.score("self", self_gff, genes, pseudo_spans, args.outdir)
    for key in MUST_BE_ONE:
        all_passed &= check(f"{key} = 100%", metrics[key] == 1.0, f"{100 * metrics[key]:.2f}%")
    for key in MUST_BE_ZERO:
        all_passed &= check(f"{key} = 0", metrics[key] == 0, str(metrics[key]))
    all_passed &= check("every reference gene matched", metrics["matched"] == len(genes),
                        f"{metrics['matched']} / {len(genes)}")
    # Alternative-start CDS share a stop with another reference CDS, so the copy has that many duplicates
    shared = sum(len(gene["cds_ids"]) - 1 for gene in genes.values())
    all_passed &= check(f"duplicate_stops = {shared} (alternative-start reference CDS)",
                        metrics["duplicate_stops"] == shared, str(metrics["duplicate_stops"]))

    # 2. Delete 10 genes
    print(f"\nCheck 2: {DELETE_COUNT} random reference genes deleted (seed {DELETE_SEED})")
    deleted_keys = random.Random(DELETE_SEED).sample(sorted(genes), DELETE_COUNT)
    drop_ids = {cds_id for key in deleted_keys for cds_id in genes[key]["cds_ids"]}
    deleted_gff = args.outdir / "reference_minus_10.gff3"
    write_reference_with_ec(args.reference_gff, args.reference_function, deleted_gff, drop_ids)

    deleted_rows = []
    for key in deleted_keys:
        gene = genes[key]
        deleted_rows.append({"seqid": key[0], "strand": key[1], "stop": key[2],
                             "locus_tag": ";".join(sorted(gene["locus_tags"])),
                             "symbols": ";".join(sorted(gene["symbols"])), "cds_ids": ";".join(gene["cds_ids"])})
        print(f"  deleted {deleted_rows[-1]['locus_tag']:<8} {deleted_rows[-1]['symbols']:<14} "
              f"{key[1]} stop {key[2]}")
    score_gff.write_tsv(args.outdir / "deleted_genes.tsv", deleted_rows)

    metrics = score_gff.score("minus10", deleted_gff, genes, pseudo_spans, args.outdir)
    with open(args.outdir / "minus10_genes.tsv") as fh:
        missed = {(row["seqid"], row["strand"], int(row["stop"]))
                  for row in csv.DictReader(fh, delimiter="\t") if row["status"] == "missed"}
    expected = set(deleted_keys)
    all_passed &= check("missed genes are exactly the deleted ones", missed == expected,
                        f"{len(missed)} missed; not deleted but missed: {len(missed - expected)}; "
                        f"deleted but not missed: {len(expected - missed)}")
    all_passed &= check(f"recall = {len(genes) - DELETE_COUNT}/{len(genes)}",
                        metrics["matched"] == len(genes) - DELETE_COUNT,
                        f"{metrics['matched']} matched, recall {100 * metrics['recall']:.2f}%")
    for key in ["precision", "start_exact", "symbol_match", "product_match", "agrees"]:
        all_passed &= check(f"{key} still 100%", metrics[key] == 1.0, f"{100 * metrics[key]:.2f}%")
    all_passed &= check("false_positives = 0", metrics["false_positives"] == 0, str(metrics["false_positives"]))

    print("\nALL CHECKS PASSED" if all_passed else "\nSOME CHECKS FAILED")
    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
