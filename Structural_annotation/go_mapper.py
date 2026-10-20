import csv
from pathlib import Path


class GOMapper:
    """Map Pfam domains found in step 2 to Gene Ontology (GO) terms."""

    COLUMNS = ["protein", "go_id", "go_name", "pfam_accession", "pfam_family"]

    def __init__(self, pfam_hits_path, pfam2go_path, output_dir, prefix="genome"):
        """
        pfam_hits_path: domain table from DomainAnnotator (..._pfam_hits.tsv)
        pfam2go_path:   pfam2go file from the Gene Ontology website
        output_dir:     folder where the results table is written
        prefix:         start of the output file name, e.g. "ecoli" -> ecoli_go_terms.tsv
        """
        self.pfam_hits_path = Path(pfam_hits_path)
        self.pfam2go_path = Path(pfam2go_path)
        self.output_dir = Path(output_dir)
        self.prefix = prefix
        self.pfam_to_go = {}
        self.domains = []
        self.go_terms = []

    def load_pfam2go(self):
        """Read the pfam2go file into {Pfam accession: [(GO id, GO name), ...]}."""
        if not self.pfam2go_path.is_file():
            raise FileNotFoundError(f"pfam2go file not found: {self.pfam2go_path}")

        self.pfam_to_go = {}

        with self.pfam2go_path.open() as fh:
            for line in fh:
                if line.startswith("!") or " > " not in line:
                    continue  # skip comments and empty lines

                left, right = line.strip().split(" > ")
                pfam_accession = left.split()[0].replace("Pfam:", "")
                go_name, go_id = right.rsplit(" ; ", 1)
                go_name = go_name.replace("GO:", "", 1)

                self.pfam_to_go.setdefault(pfam_accession, []).append((go_id, go_name))

        print(f"Loaded GO terms for {len(self.pfam_to_go)} Pfam families")

    def load_domains(self):
        """Read the step 2 domain table."""
        if not self.pfam_hits_path.is_file():
            raise FileNotFoundError(f"Pfam hits file not found: {self.pfam_hits_path}")

        with self.pfam_hits_path.open() as fh:
            self.domains = list(csv.DictReader(fh, delimiter="\t"))

        print(f"Loaded {len(self.domains)} domains from {self.pfam_hits_path.name}")

    def map_terms(self):
        """Give each protein the GO terms of its Pfam domains (each term listed once per protein)."""
        self.go_terms = []
        seen = set()

        for domain in self.domains:
            # Step 2 keeps the version ("PF00696.34"); pfam2go does not ("PF00696")
            pfam_accession = domain["pfam_accession"].split(".")[0]

            for go_id, go_name in self.pfam_to_go.get(pfam_accession, []):
                key = (domain["protein"], go_id)
                if key in seen:
                    continue
                seen.add(key)

                self.go_terms.append({
                    "protein": domain["protein"],
                    "go_id": go_id,
                    "go_name": go_name,
                    "pfam_accession": pfam_accession,
                    "pfam_family": domain["pfam_family"],
                })

        print(f"Mapped {len(self.go_terms)} GO terms")

    def write_output(self):
        """Write one row per protein–GO term pair to a tab-separated file."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        output_file = self.output_dir / f"{self.prefix}_go_terms.tsv"

        with output_file.open("w") as out:
            out.write("\t".join(self.COLUMNS) + "\n")
            for row in self.go_terms:
                out.write("\t".join(row[column] for column in self.COLUMNS) + "\n")

        print(f"Wrote {output_file}")
        return output_file

    def proteins_with_go(self):
        """Number of proteins that received at least one GO term."""
        return len({row["protein"] for row in self.go_terms})

    def run(self):
        """Run all steps in order and return the output file path."""
        self.load_pfam2go()
        self.load_domains()
        self.map_terms()
        print(f"Proteins with at least one GO term: {self.proteins_with_go()}")
        return self.write_output()
