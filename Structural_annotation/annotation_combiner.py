import csv
from pathlib import Path


class AnnotationCombiner:
    """Combine gene positions, homology names, Pfam domains and GO terms into one table."""

    COLUMNS = [
        "protein", "contig", "start", "end", "strand", "length_aa",
        "product", "product_source", "homology_status", "homology_identity", "homology_match",
        "pfam_families", "pfam_accessions", "go_ids", "go_names",
    ]

    def __init__(self, gff_path, products_path, pfam_hits_path, go_terms_path, output_dir, prefix="genome"):
        """
        gff_path:       gene coordinates from GenePredictor (the .gff file)
        products_path:  protein names from HomologyAnnotator (..._products.tsv)
        pfam_hits_path: domains from DomainAnnotator (..._pfam_hits.tsv)
        go_terms_path:  GO terms from GOMapper (..._go_terms.tsv)
        output_dir:     folder where the combined table is written
        prefix:         start of the output file name, e.g. "ecoli_swissprot" -> ecoli_swissprot_annotation.tsv
        """
        self.gff_path = Path(gff_path)
        self.products_path = Path(products_path)
        self.pfam_hits_path = Path(pfam_hits_path)
        self.go_terms_path = Path(go_terms_path)
        self.output_dir = Path(output_dir)
        self.prefix = prefix
        self.genes = {}
        self.products = {}
        self.domains = {}
        self.go_terms = {}
        self.annotations = []

    @staticmethod
    def read_tsv(path):
        """Read a tab-separated file with a header row into a list of dicts."""
        path = Path(path)
        if not path.is_file():
            raise FileNotFoundError(f"File not found: {path}")
        with path.open() as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    def load_genes(self):
        """Read gene positions from the GFF, keyed by protein ID."""
        self.genes = {}

        with self.gff_path.open() as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                contig, _, feature, start, end, _, strand, _, attributes = line.rstrip("\n").split("\t")
                if feature != "CDS":
                    continue
                protein_id = attributes.split("ID=")[1].split(";")[0]
                self.genes[protein_id] = {
                    "contig": contig,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "length_aa": (int(end) - int(start) + 1) // 3 - 1,  # minus the stop codon (approximate for partial genes)
                }

        print(f"Loaded {len(self.genes)} genes")

    def load_evidence(self):
        """Read the homology, Pfam and GO results, grouped by protein."""
        self.products = {row["protein"]: row for row in self.read_tsv(self.products_path)}

        self.domains = {}
        for row in self.read_tsv(self.pfam_hits_path):
            self.domains.setdefault(row["protein"], []).append(row)
        for rows in self.domains.values():
            rows.sort(key=lambda row: int(row["env_from"]))  # order domains along the protein

        self.go_terms = {}
        for row in self.read_tsv(self.go_terms_path):
            self.go_terms.setdefault(row["protein"], []).append(row)

        print(f"Loaded names for {len(self.products)}, domains for {len(self.domains)}, "
              f"GO terms for {len(self.go_terms)} proteins")

    @staticmethod
    def unique(values):
        """Remove repeats but keep the original order."""
        return list(dict.fromkeys(values))

    @staticmethod
    def name_from_domain(domains):
        """Describe a protein by its strongest domain, e.g. 'ABC transporter domain-containing protein'."""
        best = max(domains, key=lambda row: float(row["domain_score"]))
        if best["pfam_family"].startswith("DUF"):
            return f"{best['pfam_family']} domain-containing protein"
        description = best["description"].removesuffix(" domain")  # avoid "... domain domain-containing"
        return f"{description} domain-containing protein"

    def choose_product(self, homology, domains):
        """Pick the protein name and where it came from, using the simplest rule that works."""
        status = homology["status"] if homology else "no_match"

        if status == "confident":
            return homology["product"], "homology"

        domain_name = self.name_from_domain(domains) if domains else None
        # A match to an "Uncharacterized ..." protein says less than a known domain, so the domain wins,
        # unless the domain is itself of unknown function (DUF, UPF): then the homology name is more specific
        domain_unknown = domain_name is None or domain_name.startswith("DUF") or any(
            words in domain_name.lower() for words in ("unknown function", "uncharacterised", "uncharacterized"))
        if status == "uncharacterized_match" and domain_unknown:
            return homology["product"], "homology"
        if domain_name:
            return domain_name, "pfam_domain"
        return "hypothetical protein", "none"

    def combine(self):
        """Build one row per gene with all the evidence side by side."""
        self.annotations = []

        for protein_id, gene in self.genes.items():
            homology = self.products.get(protein_id)
            domains = self.domains.get(protein_id, [])
            go_terms = self.go_terms.get(protein_id, [])
            product, source = self.choose_product(homology, domains)

            self.annotations.append({
                "protein": protein_id,
                **gene,
                "product": product,
                "product_source": source,
                "homology_status": homology["status"] if homology else "no_match",
                "homology_identity": homology["identity"] if homology else "",
                "homology_match": homology["match_id"] if homology else "",
                "pfam_families": ";".join(self.unique(row["pfam_family"] for row in domains)),
                "pfam_accessions": ";".join(self.unique(row["pfam_accession"] for row in domains)),
                "go_ids": ";".join(self.unique(row["go_id"] for row in go_terms)),
                "go_names": ";".join(self.unique(row["go_name"] for row in go_terms)),
            })

    def write_output(self):
        """Write the combined table, one row per gene."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        output_file = self.output_dir / f"{self.prefix}_annotation.tsv"

        with output_file.open("w") as out:
            out.write("\t".join(self.COLUMNS) + "\n")
            for row in self.annotations:
                out.write("\t".join(str(row[column]) for column in self.COLUMNS) + "\n")

        print(f"Wrote {output_file}")
        return output_file

    def source_counts(self):
        """How many genes got their name from each source."""
        counts = {}
        for row in self.annotations:
            counts[row["product_source"]] = counts.get(row["product_source"], 0) + 1
        return counts

    def run(self):
        """Run all steps in order and return the output file path."""
        self.load_genes()
        self.load_evidence()
        self.combine()
        print(f"Product name sources: {self.source_counts()}")
        return self.write_output()
