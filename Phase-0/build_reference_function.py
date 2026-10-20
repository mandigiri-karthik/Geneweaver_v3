"""Join the frozen RefSeq reference CDS to UniProt, to add EC numbers the GFF does not carry.

Usage (from Phase-0/):
    python3 build_reference_function.py \
        --gff reference/ncbi_dataset/data/GCF_000005845.2/genomic.gff \
        --uniprot reference/uniprot/UP000000625_ecoli_k12.tsv \
        --output reference/reference_function.tsv

One row per reference CDS (unique GFF ID). Each CDS is linked to a UniProt entry by:
    1. the UniProtKB/Swiss-Prot accession in its Dbxref, or
    2. its b-number locus tag, found in UniProt's ordered locus names (fallback)
EC numbers, the UniProt gene name and protein name come from that entry; the RefSeq
gene symbol and product are kept alongside so the scorer can use either.

Uses only the Python standard library.
"""
import argparse
import csv
import re
from pathlib import Path

COLUMNS = [
    "cds_id", "locus_tag", "strand", "start", "end", "pseudo",
    "refseq_gene", "refseq_product", "uniprot_accession", "link_method",
    "uniprot_gene", "uniprot_protein_name", "ec_numbers",
]


def parse_attributes(text):
    return dict(item.split("=", 1) for item in text.strip(";").split(";") if "=" in item)


def read_reference_cds(gff_path):
    """Return {CDS ID: record}; split CDS rows are merged into one span."""
    cds = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attributes = parse_attributes(fields[8])
            start, end = int(fields[3]), int(fields[4])
            record = cds.get(attributes["ID"])
            if record:
                record["start"] = min(record["start"], start)
                record["end"] = max(record["end"], end)
                continue
            swissprot = re.search(r"UniProtKB/Swiss-Prot:([A-Z0-9]+)", attributes.get("Dbxref", ""))
            cds[attributes["ID"]] = {
                "cds_id": attributes["ID"],
                "locus_tag": attributes.get("locus_tag", ""),
                "strand": fields[6],
                "start": start,
                "end": end,
                "pseudo": "true" if attributes.get("pseudo") == "true" else "false",
                "refseq_gene": attributes.get("gene", ""),
                "refseq_product": attributes.get("product", ""),
                "dbxref_accession": swissprot.group(1) if swissprot else "",
            }
    return cds


def read_uniprot(tsv_path):
    """Return ({accession: entry}, {b-number: [accession, ...]})."""
    by_accession, by_locus = {}, {}
    with open(tsv_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            by_accession[row["Entry"]] = row
            for name in row["Gene Names (ordered locus)"].split():
                if re.fullmatch(r"b\d{4}", name):
                    by_locus.setdefault(name, []).append(row["Entry"])
    return by_accession, by_locus


def build(gff_path, uniprot_path, output_path):
    cds = read_reference_cds(gff_path)
    by_accession, by_locus = read_uniprot(uniprot_path)

    counts = {"dbxref": 0, "locus_tag": 0, "none": 0, "ambiguous_locus": 0}
    rows = []
    for record in cds.values():
        accession, method = record["dbxref_accession"], "dbxref"
        if accession not in by_accession:
            candidates = by_locus.get(record["locus_tag"], [])
            if len(candidates) == 1:
                accession, method = candidates[0], "locus_tag"
            else:
                if len(candidates) > 1:
                    counts["ambiguous_locus"] += 1
                accession, method = "", "none"
        counts[method] += 1

        entry = by_accession.get(accession, {})
        ec_numbers = [ec.strip() for ec in entry.get("EC number", "").split(";") if ec.strip()]
        rows.append({
            **{key: record[key] for key in COLUMNS if key in record},
            "uniprot_accession": accession,
            "link_method": method,
            "uniprot_gene": entry.get("Gene Names (primary)", ""),
            "uniprot_protein_name": entry.get("Protein names", ""),
            "ec_numbers": ";".join(ec_numbers),
        })

    output_path = Path(output_path)
    with output_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    with_ec = sum(1 for row in rows if row["ec_numbers"])
    print(f"Wrote {len(rows)} reference CDS to {output_path}")
    print(f"Linked by Dbxref: {counts['dbxref']}, by locus tag: {counts['locus_tag']}, "
          f"unlinked: {counts['none']} (ambiguous locus tag: {counts['ambiguous_locus']})")
    print(f"CDS with at least one EC number: {with_ec}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--uniprot", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    build(args.gff, args.uniprot, args.output)


if __name__ == "__main__":
    main()
