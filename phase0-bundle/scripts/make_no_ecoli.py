#!/usr/bin/env python3
"""Derive the leakage-control Swiss-Prot FASTA by removing Escherichia and Shigella.

Reads a gzipped Swiss-Prot FASTA and writes a gzipped copy with every entry whose
`OS=` organism field names Escherichia or Shigella removed. Standard library only.

Phase 0 reference counts: 575,748 input -> 545,523 output (30,225 removed).
"""
from __future__ import annotations

import argparse
import gzip
import re
import sys

# OS= runs until the next KEY= token or end of line.
_OS = re.compile(r"\bOS=(.*?)(?=\s+[A-Z]{2}=|$)")
_EXCLUDE = re.compile(r"\b(Escherichia|Shigella)\b", re.IGNORECASE)


def excluded(header: str) -> bool:
    m = _OS.search(header)
    return bool(m) and bool(_EXCLUDE.search(m.group(1)))


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="uniprot_sprot.fasta.gz")
    p.add_argument("--output", required=True, help="uniprot_sprot_no_ecoli.fasta.gz")
    p.add_argument("--expect-kept", type=int, default=None,
                   help="fail if the kept count differs (Phase 0 used 545523)")
    a = p.parse_args()

    kept = removed = 0
    keep_current = True
    # mtime=0 so the output is byte-reproducible across runs.
    with gzip.open(a.input, "rt") as fin, gzip.open(a.output, "wt", compresslevel=9, mtime=0) as fout:
        for line in fin:
            if line.startswith(">"):
                keep_current = not excluded(line)
                kept += keep_current
                removed += not keep_current
            if keep_current:
                fout.write(line)

    print(f"kept {kept}  removed {removed}  total {kept + removed}")
    if a.expect_kept is not None and kept != a.expect_kept:
        print(f"ERROR: expected {a.expect_kept} kept, got {kept}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
