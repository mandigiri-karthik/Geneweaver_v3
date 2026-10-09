"""Download a reference genome (FASTA) and its RefSeq annotation (GFF3) from NCBI.

Default: Lactobacillus sakei subsp. sakei 23K, chromosome NC_007576.1 - the
genome AGMIAL was originally used to annotate (Chaillou et al., Nat Biotechnol 2005).

    python scripts/download_data.py                      # L. sakei 23K
    python scripts/download_data.py --accession NC_000913.3 --name ecoli_k12

Run this on your own machine (needs internet access to eutils.ncbi.nlm.nih.gov).
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import requests

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
GFF_URL = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi"


def fetch(url: str, params: dict, dest: Path) -> None:
    for attempt in range(3):
        r = requests.get(url, params=params, timeout=120)
        if r.ok and r.text.strip():
            dest.write_text(r.text)
            print(f"  wrote {dest} ({len(r.text):,} bytes)")
            return
        time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"Download failed: {r.url} (HTTP {r.status_code})")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession", default="NC_007576.1")
    ap.add_argument("--name", default="lsakei_23k")
    ap.add_argument("--outdir", default="data/raw")
    a = ap.parse_args(argv)

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {a.accession} -> {out}/{a.name}.*")
    fetch(EUTILS, {"db": "nuccore", "id": a.accession, "rettype": "fasta", "retmode": "text"},
          out / f"{a.name}.fna")
    fetch(GFF_URL, {"db": "nuccore", "report": "gff3", "id": a.accession, "retmode": "text"},
          out / f"{a.name}.gff3")


if __name__ == "__main__":
    main()
