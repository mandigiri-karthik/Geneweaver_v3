from pathlib import Path

import pyhmmer


class DomainAnnotator:
    """Search predicted proteins against an HMM database (e.g. Pfam) with HMMER."""

    COLUMNS = [
        "protein", "pfam_family", "pfam_accession", "description",
        "env_from", "env_to", "domain_evalue", "domain_score",
    ]

    def __init__(self, protein_path, hmm_db_path, output_dir, prefix="genome",
                 use_gathering=True, evalue=1e-5, cpus=0):
        """
        protein_path:  protein FASTA from GenePredictor (the .faa file)
        hmm_db_path:   HMM database, e.g. Pfam-A.hmm
        output_dir:    folder where the results table is written
        prefix:        start of the output file name, e.g. "ecoli" -> ecoli_pfam_hits.tsv
        use_gathering: True to use the database's curated per-family cutoffs (Pfam has them)
        evalue:        E-value cutoff, used only when use_gathering is False
        cpus:          number of CPU cores; 0 means use all
        """
        self.protein_path = Path(protein_path)
        self.hmm_db_path = Path(hmm_db_path)
        self.output_dir = Path(output_dir)
        self.prefix = prefix
        self.use_gathering = use_gathering
        self.evalue = evalue
        self.cpus = cpus
        self.alphabet = pyhmmer.easel.Alphabet.amino()
        self.proteins = None
        self.domains = []

    def load_proteins(self):
        """Read all protein sequences into memory."""
        if not self.protein_path.is_file():
            raise FileNotFoundError(f"Protein FASTA not found: {self.protein_path}")

        with pyhmmer.easel.SequenceFile(self.protein_path, digital=True, alphabet=self.alphabet) as seq_file:
            self.proteins = seq_file.read_block()

        print(f"Loaded {len(self.proteins)} proteins from {self.protein_path.name}")

    def search(self):
        """Search every HMM in the database against the proteins and keep significant domains."""
        if not self.hmm_db_path.is_file():
            raise FileNotFoundError(f"HMM database not found: {self.hmm_db_path}")

        if self.use_gathering:
            options = {"bit_cutoffs": "gathering"}
        else:
            options = {"E": self.evalue, "domE": self.evalue}

        self.domains = []

        with pyhmmer.plan7.HMMFile(self.hmm_db_path) as hmm_file:
            for hits in pyhmmer.hmmsearch(hmm_file, self.proteins, cpus=self.cpus, **options):
                for hit in hits.included:
                    for domain in hit.domains.included:
                        self.domains.append({
                            "protein": hit.name,
                            "pfam_family": hits.query.name,
                            "pfam_accession": hits.query.accession,
                            "description": hits.query.description,
                            "env_from": domain.env_from,
                            "env_to": domain.env_to,
                            "domain_evalue": f"{domain.i_evalue:.2e}",
                            "domain_score": f"{domain.score:.1f}",
                        })

        print(f"Found {len(self.domains)} domains")

    def write_output(self):
        """Write one row per domain to a tab-separated file."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        output_file = self.output_dir / f"{self.prefix}_pfam_hits.tsv"

        with output_file.open("w") as out:
            out.write("\t".join(self.COLUMNS) + "\n")
            for row in self.domains:
                out.write("\t".join(str(row[column]) for column in self.COLUMNS) + "\n")

        print(f"Wrote {output_file}")
        return output_file

    def proteins_with_domains(self):
        """Number of proteins that have at least one domain."""
        return len({row["protein"] for row in self.domains})

    def run(self):
        """Run all steps in order and return the output file path."""
        self.load_proteins()
        self.search()
        print(f"Proteins with at least one domain: {self.proteins_with_domains()} / {len(self.proteins)}")
        return self.write_output()
