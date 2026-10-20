from pathlib import Path

import pyrodigal
from Bio import SeqIO


class GenePredictor:
    """Predict protein-coding genes in a genome with Pyrodigal and write the results."""

    def __init__(self, genome_path, output_dir, prefix="genome", meta=False, closed=False):
        """
        genome_path: FASTA file with one or more contigs
        output_dir:  folder where the output files are written
        prefix:      start of each output file name, e.g. "ecoli" -> ecoli_proteins.faa
        meta:        True for metagenomes or genomes under ~20 kb (no training step)
        closed:      True to stop genes running off the ends of contigs
        """
        self.genome_path = Path(genome_path)
        self.output_dir = Path(output_dir)
        self.prefix = prefix
        self.finder = pyrodigal.GeneFinder(meta=meta, closed=closed)
        self.records = []
        self.predictions = []

    def load_genome(self):
        """Read every contig from the FASTA file."""
        if not self.genome_path.is_file():
            raise FileNotFoundError(f"Genome FASTA not found: {self.genome_path}")

        self.records = list(SeqIO.parse(self.genome_path, "fasta"))

        if not self.records:
            raise ValueError(f"No sequences found in: {self.genome_path}")

        print(f"Loaded {len(self.records)} contig(s) from {self.genome_path.name}")

    def train(self):
        """Learn genome-specific gene features. Skipped in meta mode."""
        if self.finder.meta:
            print("Meta mode: using pre-trained models, no training needed")
            return

        sequences = [bytes(record.seq) for record in self.records]
        self.finder.train(*sequences)
        print("Training finished")

    def predict(self):
        """Find genes on each contig."""
        self.predictions = []

        for record in self.records:
            genes = self.finder.find_genes(bytes(record.seq))
            self.predictions.append((record.id, genes))
            print(f"{record.id}: {len(genes)} genes")

    def write_outputs(self):
        """Write proteins (.faa), gene sequences (.fna) and coordinates (.gff)."""
        self.output_dir.mkdir(parents=True, exist_ok=True)

        protein_file = self.output_dir / f"{self.prefix}_proteins.faa"
        gene_file = self.output_dir / f"{self.prefix}_genes.fna"
        gff_file = self.output_dir / f"{self.prefix}_genes.gff"

        with (
            protein_file.open("w") as proteins,
            gene_file.open("w") as nucleotides,
            gff_file.open("w") as gff,
        ):
            gff.write("##gff-version 3\n")

            for contig_id, genes in self.predictions:
                genes.write_translations(proteins, sequence_id=contig_id, include_stop=False)
                genes.write_genes(nucleotides, sequence_id=contig_id)
                genes.write_gff(gff, sequence_id=contig_id, header=False)

        print(f"Wrote {protein_file.name}, {gene_file.name}, {gff_file.name} to {self.output_dir}")
        return {"proteins": protein_file, "genes": gene_file, "gff": gff_file}

    def total_genes(self):
        """Number of genes predicted across all contigs."""
        return sum(len(genes) for _, genes in self.predictions)

    def run(self):
        """Run all steps in order and return the output file paths."""
        self.load_genome()
        self.train()
        self.predict()
        print(f"Total genes: {self.total_genes()}")
        return self.write_outputs()
