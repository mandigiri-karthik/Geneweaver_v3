from pathlib import Path

from Bio import SeqIO
import pyrodigal


def main():
    # Directory containing this script: Structural_annotation
    base_dir = Path(__file__).resolve().parent

    # Input genome is in the sibling directory: E-coli_data
    input_file = (
        base_dir.parent
        / "E-coli_data"
        / "ncbi_dataset"
        / "data"
        / "GCF_000005845.2"
        / "GCF_000005845.2_ASM584v2_genomic.fna"
    )

    # Save outputs inside Structural_annotation/pyrodigal_results.
    output_dir = base_dir / "pyrodigal_results"

    if not input_file.is_file():
        raise FileNotFoundError(f"Input FASTA not found: {input_file}")

    print(f"Reading genome: {input_file}")

    records = list(SeqIO.parse(input_file, "fasta"))

    if not records:
        raise ValueError(f"No FASTA records found in: {input_file}")

    sequences = [bytes(record.seq) for record in records]

    # Train a genome-specific model using the input sequences.
    print("Training Pyrodigal...")
    finder = pyrodigal.GeneFinder(meta=False)
    finder.train(*sequences)

    output_dir.mkdir(parents=True, exist_ok=True)

    protein_file = output_dir / "ecoli_proteins.faa"
    gene_file = output_dir / "ecoli_genes.fna"
    gff_file = output_dir / "ecoli_genes.gff"

    total_genes = 0

    print("Predicting genes...")

    with (
        protein_file.open("w") as proteins,
        gene_file.open("w") as nucleotides,
        gff_file.open("w") as gff,
    ):
        gff.write("##gff-version 3\n")

        for record, sequence in zip(records, sequences):
            genes = finder.find_genes(sequence)
            total_genes += len(genes)

            genes.write_translations(
                proteins,
                sequence_id=record.id,
                include_stop=False,
                full_id=True,
            )

            genes.write_genes(
                nucleotides,
                sequence_id=record.id,
                full_id=True,
            )

            genes.write_gff(
                gff,
                sequence_id=record.id,
                header=False,
                full_id=True,
            )

            print(f"{record.id}: {len(genes)} predicted genes")

    print(f"\nTotal predicted genes: {total_genes}")
    print(f"Proteins: {protein_file}")
    print(f"Gene DNA: {gene_file}")
    print(f"Gene coordinates: {gff_file}")


if __name__ == "__main__":
    main()