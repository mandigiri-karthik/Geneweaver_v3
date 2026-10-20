import argparse
from pathlib import Path

from Structural_annotation.annotation_combiner import AnnotationCombiner
from Structural_annotation.go_mapper import GOMapper
from Structural_annotation.homology_annotator import HomologyAnnotator
from Structural_annotation.pfam_annotator import DomainAnnotator
from Structural_annotation.predict_genes import GenePredictor

ROOT = Path(__file__).resolve().parent


def parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate a genome: predict genes, name proteins by homology, "
                    "find Pfam domains, map them to GO terms and combine the results."
    )
    parser.add_argument("genome", type=Path, help="Path to the genome FASTA file")
    parser.add_argument("--prefix", default="genome", help="Prefix for the output file names")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "Structural_annotation/pyrodigal_results",
        help="Folder for the gene prediction files",
    )
    parser.add_argument(
        "--diamond-db",
        type=Path,
        default=ROOT / "Structural_annotation/databases/swissprot",
        help="DIAMOND database for the homology search (without .dmnd)",
    )
    parser.add_argument(
        "--homology-output-dir",
        type=Path,
        default=ROOT / "Structural_annotation/homology_results",
        help="Folder for the homology results",
    )
    parser.add_argument(
        "--hmm-db",
        type=Path,
        default=ROOT / "Structural_annotation/databases/Pfam-A.hmm",
        help="HMM database for the domain search",
    )
    parser.add_argument(
        "--pfam-clans",
        type=Path,
        default=ROOT / "Structural_annotation/databases/Pfam-A.clans.tsv.gz",
        help="Pfam clan table, used to remove overlapping domains from the same clan",
    )
    parser.add_argument(
        "--domain-output-dir",
        type=Path,
        default=ROOT / "Structural_annotation/pfam_results",
        help="Folder for the domain search results",
    )
    parser.add_argument(
        "--pfam2go",
        type=Path,
        default=ROOT / "Structural_annotation/databases/pfam2go",
        help="Pfam to GO mapping file",
    )
    parser.add_argument(
        "--go-output-dir",
        type=Path,
        default=ROOT / "Structural_annotation/go_results",
        help="Folder for the GO mapping results",
    )
    parser.add_argument(
        "--annotation-output-dir",
        type=Path,
        default=ROOT / "Structural_annotation/annotation_results",
        help="Folder for the combined annotation table",
    )
    parser.add_argument(
        "--genes-only",
        action="store_true",
        help="Run gene prediction only",
    )
    parser.add_argument(
        "--reuse-domains",
        action="store_true",
        help="Reuse existing Pfam and GO results instead of recomputing them "
             "(only valid if the genome and gene prediction have not changed)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Step 1: gene prediction
    predictor = GenePredictor(
        genome_path=args.genome,
        output_dir=args.output_dir,
        prefix=args.prefix,
    )
    outputs = predictor.run()

    if args.genes_only:
        return outputs

    # Step 2: protein function by homology
    homology = HomologyAnnotator(
        protein_path=outputs["proteins"],
        diamond_db_path=args.diamond_db,
        output_dir=args.homology_output_dir,
        prefix=args.prefix,
    )
    outputs["products"] = homology.run()

    # Steps 3 and 4 do not depend on the homology database, so they can be reused
    pfam_file = args.domain_output_dir / f"{args.prefix}_pfam_hits.tsv"
    go_file = args.go_output_dir / f"{args.prefix}_go_terms.tsv"

    if args.reuse_domains and pfam_file.is_file() and go_file.is_file():
        print(f"Reusing {pfam_file.name} and {go_file.name}")
        outputs["pfam_hits"] = pfam_file
        outputs["go_terms"] = go_file
    else:
        # Step 3: domain search on the proteins from step 1
        annotator = DomainAnnotator(
            protein_path=outputs["proteins"],
            hmm_db_path=args.hmm_db,
            output_dir=args.domain_output_dir,
            prefix=args.prefix,
            clans_path=args.pfam_clans,
        )
        outputs["pfam_hits"] = annotator.run()

        # Step 4: Pfam domains -> GO terms
        go_mapper = GOMapper(
            pfam_hits_path=outputs["pfam_hits"],
            pfam2go_path=args.pfam2go,
            output_dir=args.go_output_dir,
            prefix=args.prefix,
        )
        outputs["go_terms"] = go_mapper.run()

    # Step 5: combine all evidence into one table
    combiner = AnnotationCombiner(
        gff_path=outputs["gff"],
        products_path=outputs["products"],
        pfam_hits_path=outputs["pfam_hits"],
        go_terms_path=outputs["go_terms"],
        output_dir=args.annotation_output_dir,
        prefix=f"{args.prefix}_{args.diamond_db.name}",
    )
    outputs["annotation"] = combiner.run()

    return outputs


if __name__ == "__main__":
    main()
