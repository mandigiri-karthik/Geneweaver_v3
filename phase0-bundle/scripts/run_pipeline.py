#!/usr/bin/env python3
"""Driver for the deterministic annotation pipeline (arm B).

The four annotator modules are library classes with no __main__ block - the
original driver was a notebook (`Structural_annotation/stage0.ipynb`, deleted in
d0789ef). This replaces it with a reproducible, parameterised entrypoint.

Produces, for each homology condition, exactly the filenames Phase 0's
export_pipeline_gff3.py expects:

    <prefix>_<db>_hits.tsv        <prefix>_<db>_products.tsv
    <prefix>_pfam_hits.tsv        <prefix>_go_terms.tsv
    <prefix>_<db>_annotation.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

DEFAULT_CONDITIONS = ["swissprot", "swissprot_no_ecoli"]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--repo", type=Path, required=True,
                   help="repository root (contains Structural_annotation/ and Phase-0/)")
    p.add_argument("--databases", type=Path, required=True,
                   help="directory holding swissprot.dmnd, Pfam-A.hmm, pfam2go")
    p.add_argument("--prefix", default="ecoli")
    p.add_argument("--conditions", nargs="+", default=DEFAULT_CONDITIONS,
                   help="DIAMOND database basenames to run")
    p.add_argument("--threads", type=int, default=8, help="DIAMOND threads")
    p.add_argument("--cpus", type=int, default=0, help="PyHMMER cpus (0 = all)")
    p.add_argument("--conda-env", default=None,
                   help="run DIAMOND via 'conda run -n ENV'; omit when diamond is on PATH")
    p.add_argument("--min-identity", type=float, default=40.0)
    p.add_argument("--min-coverage", type=float, default=80.0)
    a = p.parse_args()

    sa = a.repo / "Structural_annotation"
    sys.path.insert(0, str(sa))
    from annotation_combiner import AnnotationCombiner      # noqa: E402
    from go_mapper import GOMapper                          # noqa: E402
    from homology_annotator import HomologyAnnotator        # noqa: E402
    from pfam_annotator import DomainAnnotator              # noqa: E402

    proteins = sa / "pyrodigal_results" / f"{a.prefix}_proteins.faa"
    genes_gff = sa / "pyrodigal_results" / f"{a.prefix}_genes.gff"
    for f in (proteins, genes_gff):
        if not f.is_file():
            print(f"missing {f}\nRun predict_genes.py first.", file=sys.stderr)
            return 1

    hom_dir = sa / "homology_results"
    pfam_dir = sa / "pfam_results"
    go_dir = sa / "go_results"
    ann_dir = sa / "annotation_results"

    # --- domain + GO evidence: condition-independent, so computed once -------
    print("== Pfam domains (PyHMMER, gathering cutoffs) ==")
    DomainAnnotator(
        protein_path=proteins,
        hmm_db_path=a.databases / "Pfam-A.hmm",
        output_dir=pfam_dir,
        prefix=a.prefix,
        cpus=a.cpus,
    ).run()

    print("== GO terms (pfam2go) ==")
    GOMapper(
        pfam_hits_path=pfam_dir / f"{a.prefix}_pfam_hits.tsv",
        pfam2go_path=a.databases / "pfam2go",
        output_dir=go_dir,
        prefix=a.prefix,
    ).run()

    # --- homology + combine: once per leakage condition ----------------------
    for cond in a.conditions:
        print(f"== homology: {cond} ==")
        HomologyAnnotator(
            protein_path=proteins,
            diamond_db_path=a.databases / cond,
            output_dir=hom_dir,
            prefix=a.prefix,
            min_identity=a.min_identity,
            min_coverage=a.min_coverage,
            threads=a.threads,
            conda_env=a.conda_env,
        ).run()

        print(f"== combine: {cond} ==")
        AnnotationCombiner(
            gff_path=genes_gff,
            products_path=hom_dir / f"{a.prefix}_{cond}_products.tsv",
            pfam_hits_path=pfam_dir / f"{a.prefix}_pfam_hits.tsv",
            go_terms_path=go_dir / f"{a.prefix}_go_terms.tsv",
            output_dir=ann_dir,
            prefix=f"{a.prefix}_{cond}",
        ).run()

    print("\nPipeline complete. Combined tables:")
    for cond in a.conditions:
        print("  ", ann_dir / f"{a.prefix}_{cond}_annotation.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
