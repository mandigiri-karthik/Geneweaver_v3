"""Export the deterministic pipeline's combined annotation table as GFF3, for Phase 0 scoring.

Usage:
    python3 export_pipeline_gff3.py \
        --annotation <..._annotation.tsv> \
        --genes-gff <pyrodigal ..._genes.gff> \
        --swissprot <uniprot_sprot.fasta.gz> \
        --genome <genome.fna> \
        --output <out.gff3>

One CDS line per predicted gene, with:
    product   the pipeline's chosen product name
    gene      Swiss-Prot gene symbol (GN=) of the matched protein, only when the product
              name was copied from that match (product_source == homology)
    Dbxref    Pfam accessions and GO terms found for the protein
No EC numbers: the pipeline does not assign them (Swiss-Prot FASTA headers do not carry EC).

Uses only the Python standard library.
"""
import argparse
import csv
import gzip
import re
from pathlib import Path

SOURCE = "GeneWeaver-baseline"

# Characters that must be percent-encoded inside GFF3 attribute values
GFF3_ESCAPES = {"%": "%25", "\t": "%09", "\n": "%0A", "\r": "%0D", ";": "%3B", "=": "%3D", "&": "%26", ",": "%2C"}


def escape(value):
    """Percent-encode a single GFF3 attribute value."""
    return "".join(GFF3_ESCAPES.get(char, char) for char in str(value))


def read_gene_symbols(swissprot_fasta):
    """Map Swiss-Prot accession -> gene symbol (GN=) from the FASTA headers."""
    symbols = {}
    opener = gzip.open if str(swissprot_fasta).endswith(".gz") else open
    with opener(swissprot_fasta, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            accession = line.split("|")[1]
            match = re.search(r" GN=(\S+)", line)
            if match:
                symbols[accession] = match.group(1)
    return symbols


def read_pyrodigal_gff(gff_path):
    """Read score, partial flag and start type for each predicted gene, keyed by protein ID."""
    genes = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "CDS":
                continue
            attributes = dict(item.split("=", 1) for item in fields[8].strip(";").split(";") if "=" in item)
            genes[attributes["ID"]] = {
                "score": fields[5],
                "partial": attributes.get("partial", ""),
                "start_type": attributes.get("start_type", ""),
            }
    return genes


def read_sequence_lengths(genome_fasta):
    """Return [(sequence ID, length)] in file order, for the ##sequence-region lines."""
    lengths = []
    with open(genome_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                lengths.append([line[1:].split()[0], 0])
            elif lengths:
                lengths[-1][1] += len(line.strip())
    return [tuple(item) for item in lengths]


def build_attributes(row, gene_info, symbols):
    """Build the ordered GFF3 attribute list for one gene."""
    accession = row["homology_match"].split("|")[1] if row["homology_match"].count("|") >= 2 else ""
    symbol = symbols.get(accession, "") if row["product_source"] == "homology" else ""

    dbxrefs = [f"Pfam:{acc}" for acc in row["pfam_accessions"].split(";") if acc]
    dbxrefs += [go for go in row["go_ids"].split(";") if go]   # GO IDs already look like GO:0005524

    attributes = [
        ("ID", row["protein"]),
        ("locus_tag", row["protein"]),
        ("product", row["product"]),
    ]
    if symbol:
        attributes.append(("gene", symbol))
    attributes += [
        ("product_source", row["product_source"]),
        ("homology_status", row["homology_status"]),
    ]
    if row["homology_match"]:
        attributes += [("homology_match", row["homology_match"]), ("homology_identity", row["homology_identity"])]
    if gene_info["partial"]:
        attributes.append(("partial", gene_info["partial"]))
    if gene_info["start_type"]:
        attributes.append(("start_type", gene_info["start_type"]))

    text = ";".join(f"{key}={escape(value)}" for key, value in attributes)
    if dbxrefs:
        text += ";Dbxref=" + ",".join(escape(ref) for ref in dbxrefs)
    return text


def export(annotation_tsv, genes_gff, swissprot_fasta, genome_fasta, output_gff3):
    symbols = read_gene_symbols(swissprot_fasta)
    genes = read_pyrodigal_gff(genes_gff)

    with open(annotation_tsv) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    missing = [row["protein"] for row in rows if row["protein"] not in genes]
    if missing:
        raise ValueError(f"{len(missing)} proteins in the table are not in the gene GFF, e.g. {missing[:3]}")

    output_gff3 = Path(output_gff3)
    output_gff3.parent.mkdir(parents=True, exist_ok=True)
    with output_gff3.open("w") as out:
        out.write("##gff-version 3\n")
        out.write(f"# Exported by export_pipeline_gff3.py from {Path(annotation_tsv).name}\n")
        for sequence_id, length in read_sequence_lengths(genome_fasta):
            out.write(f"##sequence-region {sequence_id} 1 {length}\n")
        for row in rows:
            gene_info = genes[row["protein"]]
            columns = [
                row["contig"], SOURCE, "CDS", row["start"], row["end"],
                gene_info["score"] or ".", row["strand"], "0",
                build_attributes(row, gene_info, symbols),
            ]
            out.write("\t".join(columns) + "\n")

    with_symbol = sum(1 for row in rows if row["product_source"] == "homology"
                      and row["homology_match"].count("|") >= 2
                      and row["homology_match"].split("|")[1] in symbols)
    print(f"Wrote {len(rows)} CDS to {output_gff3} ({with_symbol} with a gene symbol)")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--annotation", required=True, type=Path, help="pipeline ..._annotation.tsv")
    parser.add_argument("--genes-gff", required=True, type=Path, help="Pyrodigal GFF from the same run")
    parser.add_argument("--swissprot", required=True, type=Path, help="Swiss-Prot FASTA used for the homology search")
    parser.add_argument("--genome", required=True, type=Path, help="genome FASTA (for sequence lengths)")
    parser.add_argument("--output", required=True, type=Path, help="output GFF3 path")
    args = parser.parse_args()
    export(args.annotation, args.genes_gff, args.swissprot, args.genome, args.output)


if __name__ == "__main__":
    main()
