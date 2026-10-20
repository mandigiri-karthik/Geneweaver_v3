import gzip
from pathlib import Path

import pyhmmer


class DomainAnnotator:
    """Search predicted proteins against an HMM database (e.g. Pfam) with HMMER."""

    COLUMNS = [
        "protein", "pfam_family", "pfam_accession", "description",
        "env_from", "env_to", "domain_evalue", "domain_score",
    ]

    def __init__(self, protein_path, hmm_db_path, output_dir, prefix="genome",
                 use_gathering=True, evalue=1e-5, cpus=0, clans_path=None):
        """
        protein_path:  protein FASTA from GenePredictor (the .faa file)
        hmm_db_path:   HMM database, e.g. Pfam-A.hmm
        output_dir:    folder where the results table is written
        prefix:        start of the output file name, e.g. "ecoli" -> ecoli_pfam_hits.tsv
        use_gathering: True to use the database's curated per-family cutoffs (Pfam has them)
        evalue:        E-value cutoff, used only when use_gathering is False
        cpus:          number of CPU cores; 0 means use all
        clans_path:    Pfam-A.clans.tsv(.gz); if given, overlapping domains from the same clan
                       are reduced to the best-scoring one (all hits are still kept in ..._pfam_hits_all.tsv)
        """
        self.protein_path = Path(protein_path)
        self.hmm_db_path = Path(hmm_db_path)
        self.output_dir = Path(output_dir)
        self.prefix = prefix
        self.use_gathering = use_gathering
        self.evalue = evalue
        self.cpus = cpus
        self.clans_path = Path(clans_path) if clans_path else None
        self.alphabet = pyhmmer.easel.Alphabet.amino()
        self.proteins = None
        self.domains = []
        self.all_domains = []
        self.clan_of = {}

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

    def load_clans(self):
        """Read the Pfam family -> clan table. Families without a clan are their own group."""
        if not self.clans_path.is_file():
            raise FileNotFoundError(f"Pfam clans file not found: {self.clans_path}")

        opener = gzip.open if self.clans_path.suffix == ".gz" else open
        self.clan_of = {}
        with opener(self.clans_path, "rt") as fh:
            for line in fh:
                family, clan = line.rstrip("\n").split("\t")[:2]
                self.clan_of[family] = clan or family

        print(f"Loaded clans for {len(self.clan_of)} Pfam families")

    def clan(self, domain):
        """Clan of a domain's family (accession without version, e.g. PF00696.34 -> PF00696)."""
        family = domain["pfam_accession"].split(".")[0]
        return self.clan_of.get(family, family)

    def resolve_clan_overlaps(self):
        """Where domains from the same clan overlap on a protein, keep only the best-scoring one."""
        self.all_domains = self.domains
        by_protein = {}
        for domain in self.domains:
            by_protein.setdefault(domain["protein"], []).append(domain)

        self.domains = []
        for domains in by_protein.values():
            kept = []
            for domain in sorted(domains, key=lambda d: float(d["domain_score"]), reverse=True):
                overlaps_same_clan = any(
                    self.clan(k) == self.clan(domain)
                    and k["env_from"] <= domain["env_to"] and domain["env_from"] <= k["env_to"]
                    for k in kept
                )
                if not overlaps_same_clan:
                    kept.append(domain)
            self.domains.extend(sorted(kept, key=lambda d: d["env_from"]))

        removed = len(self.all_domains) - len(self.domains)
        print(f"Removed {removed} overlapping same-clan domains; {len(self.domains)} domains kept")

    def write_table(self, rows, output_file):
        """Write domains to a tab-separated file, one row per domain."""
        with output_file.open("w") as out:
            out.write("\t".join(self.COLUMNS) + "\n")
            for row in rows:
                out.write("\t".join(str(row[column]) for column in self.COLUMNS) + "\n")
        print(f"Wrote {output_file}")

    def write_output(self):
        """Write the final domains, plus every hit before overlap resolution when it was applied."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        output_file = self.output_dir / f"{self.prefix}_pfam_hits.tsv"
        self.write_table(self.domains, output_file)

        if self.all_domains:
            self.write_table(self.all_domains, self.output_dir / f"{self.prefix}_pfam_hits_all.tsv")

        return output_file

    def proteins_with_domains(self):
        """Number of proteins that have at least one domain."""
        return len({row["protein"] for row in self.domains})

    def run(self):
        """Run all steps in order and return the output file path."""
        self.load_proteins()
        self.search()
        if self.clans_path:
            self.load_clans()
            self.resolve_clan_overlaps()
        print(f"Proteins with at least one domain: {self.proteins_with_domains()} / {len(self.proteins)}")
        return self.write_output()
