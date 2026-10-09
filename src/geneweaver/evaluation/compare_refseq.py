"""Compare CAM gene calls with a reference (e.g. RefSeq/PGAP) GFF3.

Matching rule (standard for prokaryotic gene-finder benchmarks, e.g. the
Prodigal paper): two CDS are the *same gene* if they share contig, strand and
stop coordinate. If the start coordinate also agrees, the call is an *exact*
match; otherwise it is a *start error* (right gene, wrong start codon).

Reference CDS marked pseudo=true are excluded and reported separately.

Usage:
    python -m geneweaver.evaluation.compare_refseq PRED.gff3 REF.gff [--out results/cam_eval.json]
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from urllib.parse import unquote


def _attrs(col9: str) -> dict[str, str]:
    out = {}
    for part in col9.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = unquote(v)
    return out


def read_cds(gff_path: str | Path) -> tuple[list[dict], int]:
    """Return (non-pseudo CDS list, number of pseudo CDS skipped).

    Multi-segment CDS (e.g. programmed frameshifts) share an ID and are merged
    into one span.
    """
    merged: dict[str, dict] = {}
    pseudo_ids = set()
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            a = _attrs(cols[8])
            key = a.get("ID") or f"{cols[0]}:{cols[3]}-{cols[4]}{cols[6]}"
            if a.get("pseudo", "").lower() == "true":
                pseudo_ids.add(key)
                continue
            start, end = int(cols[3]), int(cols[4])
            rec = merged.setdefault(
                key,
                {"id": key, "contig": cols[0], "strand": cols[6], "start": start, "end": end,
                 "locus_tag": a.get("locus_tag", ""), "product": a.get("product", "")},
            )
            rec["start"], rec["end"] = min(rec["start"], start), max(rec["end"], end)
    return list(merged.values()), len(pseudo_ids)


def _stop(c: dict) -> int:
    return c["end"] if c["strand"] == "+" else c["start"]


def _start(c: dict) -> int:
    return c["start"] if c["strand"] == "+" else c["end"]


def compare(pred_gff: str | Path, ref_gff: str | Path) -> dict:
    pred, _ = read_cds(pred_gff)
    ref, n_pseudo = read_cds(ref_gff)

    ref_by_stop = {(c["contig"], c["strand"], _stop(c)): c for c in ref}
    pred_by_stop = {(c["contig"], c["strand"], _stop(c)): c for c in pred}

    exact, start_errors, shifts = 0, [], Counter()
    for key, p in pred_by_stop.items():
        r = ref_by_stop.get(key)
        if r is None:
            continue
        if _start(p) == _start(r):
            exact += 1
        else:
            # positive shift = CAM gene is LONGER than the reference (upstream start)
            shift = abs(_stop(p) - _start(p)) - abs(_stop(r) - _start(r))
            shifts[shift] += 1
            start_errors.append({"pred": p["id"], "ref_locus_tag": r["locus_tag"],
                                 "ref_product": r["product"], "shift_bp": shift})

    matched = exact + len(start_errors)
    missed = [r for k, r in ref_by_stop.items() if k not in pred_by_stop]
    extra = [p for k, p in pred_by_stop.items() if k not in ref_by_stop]

    n_ref, n_pred = len(ref_by_stop), len(pred_by_stop)
    return {
        "n_reference_cds": n_ref,
        "n_reference_pseudo_skipped": n_pseudo,
        "n_predicted_cds": n_pred,
        "matched_by_stop": matched,
        "exact_match": exact,
        "start_errors": len(start_errors),
        "missed": len(missed),
        "extra": len(extra),
        "sensitivity_pct": round(100 * matched / n_ref, 2) if n_ref else 0.0,
        "precision_pct": round(100 * matched / n_pred, 2) if n_pred else 0.0,
        "exact_pct_of_matched": round(100 * exact / matched, 2) if matched else 0.0,
        "start_error_pct_of_matched": round(100 * len(start_errors) / matched, 2) if matched else 0.0,
        "start_shift_longer_vs_shorter": {
            "cam_longer": sum(v for k, v in shifts.items() if k > 0),
            "cam_shorter": sum(v for k, v in shifts.items() if k < 0),
        },
        "examples": {
            "start_errors": start_errors[:20],
            "missed": [{"locus_tag": r["locus_tag"], "product": r["product"],
                        "length_bp": r["end"] - r["start"] + 1} for r in missed[:20]],
            "extra": [{"id": p["id"], "length_bp": p["end"] - p["start"] + 1} for p in extra[:20]],
        },
    }


def format_table(m: dict) -> str:
    rows = [
        ("Reference CDS (non-pseudo)", m["n_reference_cds"]),
        ("Reference pseudogene CDS skipped", m["n_reference_pseudo_skipped"]),
        ("CAM predicted CDS", m["n_predicted_cds"]),
        ("Same gene (stop match)", f"{m['matched_by_stop']}  (sensitivity {m['sensitivity_pct']}%)"),
        ("  exact (start + stop)", f"{m['exact_match']}  ({m['exact_pct_of_matched']}% of matched)"),
        ("  wrong start codon", f"{m['start_errors']}  ({m['start_error_pct_of_matched']}% of matched)"),
        ("Missed reference genes", m["missed"]),
        ("Extra CAM genes", f"{m['extra']}  (precision {m['precision_pct']}%)"),
    ]
    w = max(len(r[0]) for r in rows)
    return "\n".join(f"{k.ljust(w)}  {v}" for k, v in rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pred_gff")
    ap.add_argument("ref_gff")
    ap.add_argument("--out", default=None, help="write full metrics as JSON")
    args = ap.parse_args(argv)
    m = compare(args.pred_gff, args.ref_gff)
    print(format_table(m))
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(m, indent=2))
        print(f"\nFull metrics written to {args.out}")


if __name__ == "__main__":
    main()
