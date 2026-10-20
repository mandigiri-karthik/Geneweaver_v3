import csv
import subprocess
import re
from pathlib import Path

from Bio import SeqIO


class HomologyAnnotator:
    """Name predicted proteins by their most similar protein in Swiss-Prot (DIAMOND search)."""

    DIAMOND_FIELDS = [
        "qseqid", "sseqid", "pident", "length", "qlen", "slen",
        "qcovhsp", "scovhsp", "evalue", "bitscore", "stitle",
    ]
    COLUMNS = [
        "protein", "product", "status", "match_id", "match_organism",
        "identity", "query_coverage", "match_coverage", "evalue", "bitscore",
    ]

    def __init__(self, protein_path, diamond_db_path, output_dir, prefix="genome",
                 min_identity=40.0, min_coverage=80.0, evalue=1e-5, max_targets=5,
                 threads=8, conda_env="diamond"):
        """
        protein_path:    protein FASTA from GenePredictor (the .faa file)
        diamond_db_path: DIAMOND database, e.g. databases/swissprot (without .dmnd)
        output_dir:      folder where the results are written
        prefix:          start of the output file names, e.g. "ecoli"
        min_identity:    minimum % identity to copy a name
        min_coverage:    minimum % of both proteins covered by the alignment
        evalue:          maximum E-value for DIAMOND to report a match
        max_targets:     number of matches kept per protein
        threads:         CPU threads for DIAMOND
        conda_env:       conda environment where DIAMOND is installed
        """
        self.protein_path = Path(protein_path)
        self.diamond_db_path = Path(diamond_db_path)
        self.output_dir = Path(output_dir)
        self.prefix = prefix
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        self.evalue = evalue
        self.max_targets = max_targets
        self.threads = threads
        self.conda_env = conda_env
        self.hits_file = self.output_dir / f"{prefix}_{self.diamond_db_path.name}_hits.tsv"
        self.hits = {}
        self.annotations = []

    def run_diamond(self):
        """Search all proteins against the database and write the raw matches."""
        if not self.protein_path.is_file():
            raise FileNotFoundError(f"Protein FASTA not found: {self.protein_path}")

        self.output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            "conda", "run", "-n", self.conda_env,
            "diamond", "blastp",
            "--query", str(self.protein_path),
            "--db", str(self.diamond_db_path),
            "--out", str(self.hits_file),
            "--outfmt", "6", *self.DIAMOND_FIELDS,
            "--more-sensitive",
            "--evalue", str(self.evalue),
            "--max-target-seqs", str(self.max_targets),
            "--threads", str(self.threads),
            "--quiet",
        ]
        print(f"Running DIAMOND against {self.diamond_db_path.name}...")
        subprocess.run(cmd, check=True)
        print(f"Wrote raw matches to {self.hits_file}")

    def load_hits(self):
        """Read the raw matches into {protein: [match, ...]}, best match first, duplicates removed."""
        self.hits = {}
        seen = set()

        with self.hits_file.open() as fh:
            for values in csv.reader(fh, delimiter="\t"):
                hit = dict(zip(self.DIAMOND_FIELDS, values))
                key = (hit["qseqid"], hit["sseqid"])
                if key in seen:
                    continue  # DIAMOND sometimes repeats a row
                seen.add(key)
                self.hits.setdefault(hit["qseqid"], []).append(hit)

        print(f"Loaded matches for {len(self.hits)} proteins")

    @staticmethod
    def parse_title(title):
        """Split a Swiss-Prot title into (protein name, organism)."""
        description = title.split(" ", 1)[1]               # drop the "sp|P00561|AK1H_ECOLI" part
        name = description.split(" OS=")[0]
        organism = description.split(" OS=")[1].split(" OX=")[0] if " OS=" in description else ""
        return name, organism

    def passes(self, hit):
        """True if the match is similar enough to copy its name."""
        return (
            float(hit["pident"]) >= self.min_identity
            and float(hit["qcovhsp"]) >= self.min_coverage
            and float(hit["scovhsp"]) >= self.min_coverage
        )

    @staticmethod
    def is_informative(name):
        """False for names that say nothing about function."""
        uninformative_starts = ("Uncharacterized", "Putative uncharacterized", "UPF0")
        gene_name_only = re.fullmatch(r"(Putative )?[Pp]rotein Y[a-z]{2}[A-Z0-9]*", name)
        return not (name.startswith(uninformative_starts) or gene_name_only)


    def assign_names(self):
        """Choose one name per protein, including proteins with no match."""
        protein_ids = [record.id for record in SeqIO.parse(self.protein_path, "fasta")]
        self.annotations = []

        for protein_id in protein_ids:
            matches = self.hits.get(protein_id, [])
            good = [hit for hit in matches if self.passes(hit)]
            informative = [hit for hit in good if self.is_informative(self.parse_title(hit["stitle"])[0])]

            if informative:
                best, status = informative[0], "confident"
            elif good:
                best, status = good[0], "uncharacterized_match"
            elif matches:
                best, status = matches[0], "weak_match"
            else:
                best, status = None, "no_match"

            if best is None:
                self.annotations.append({
                    "protein": protein_id, "product": "hypothetical protein", "status": status,
                    "match_id": "", "match_organism": "", "identity": "", "query_coverage": "",
                    "match_coverage": "", "evalue": "", "bitscore": "",
                })
                continue

            name, organism = self.parse_title(best["stitle"])
            self.annotations.append({
                "protein": protein_id,
                "product": name if status != "weak_match" else "hypothetical protein",
                "status": status,
                "match_id": best["sseqid"],
                "match_organism": organism,
                "identity": best["pident"],
                "query_coverage": best["qcovhsp"],
                "match_coverage": best["scovhsp"],
                "evalue": best["evalue"],
                "bitscore": best["bitscore"],
            })

    def write_output(self):
        """Write one row per protein to a tab-separated file."""
        output_file = self.output_dir / f"{self.prefix}_{self.diamond_db_path.name}_products.tsv"

        with output_file.open("w") as out:
            out.write("\t".join(self.COLUMNS) + "\n")
            for row in self.annotations:
                out.write("\t".join(str(row[column]) for column in self.COLUMNS) + "\n")

        print(f"Wrote {output_file}")
        return output_file

    def status_counts(self):
        """How many proteins ended up in each status."""
        counts = {}
        for row in self.annotations:
            counts[row["status"]] = counts.get(row["status"], 0) + 1
        return counts

    def run(self):
        """Run all steps in order and return the output file path."""
        self.run_diamond()
        self.load_hits()
        self.assign_names()
        print(f"Status counts: {self.status_counts()}")
        return self.write_output()
