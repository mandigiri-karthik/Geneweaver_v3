#!/usr/bin/env bash
# Phase 0 end to end: predict -> annotate -> export GFF3 -> score -> sanity checks.
# Usage:  scripts/run_phase0.sh [REPO_ROOT]
set -euo pipefail

ROOT="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
DB="${GENEWEAVER_DB:-/data}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="$ROOT/Structural_annotation"
P="$ROOT/Phase-0"
GENOME="$P/reference/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"

echo "repo      : $ROOT"
echo "databases : $DB"
diamond --version
python3 -c "import pyrodigal,pyhmmer,Bio;print('pyrodigal',pyrodigal.__version__,'| pyhmmer',pyhmmer.__version__,'| biopython',Bio.__version__)"

for f in "$S/homology_annotator.py" "$S/pfam_annotator.py" "$S/go_mapper.py" "$S/annotation_combiner.py"; do
    [ -f "$f" ] || { echo "ERROR: $f missing - merge branch 'dev' into 'main' first." >&2; exit 1; }
done

step() { echo; echo "### $*"; }

step "1/5  gene prediction (Pyrodigal, translation table 11)"
python3 "$S/predict_genes.py"

step "2/5  homology, domains, GO, combine"
python3 "$HERE/run_pipeline.py" --repo "$ROOT" --databases "$DB"

step "3/5  GFF3 export"
cd "$P"
for cond in swissprot swissprot_no_ecoli; do
    db="uniprot_sprot.fasta.gz"
    [ "$cond" = "swissprot_no_ecoli" ] && db="uniprot_sprot_no_ecoli.fasta.gz"
    python3 export_pipeline_gff3.py \
        --annotation "$S/annotation_results/ecoli_${cond}_annotation.tsv" \
        --genes-gff  "$S/pyrodigal_results/ecoli_genes.gff" \
        --swissprot  "$DB/$db" \
        --genome     "$GENOME" \
        --output     "pipeline/pipeline_${cond}.gff3"
done

step "4/5  scoring against the frozen reference"
python3 score_gff.py \
    --pred bakta=bakta/result.gff3 \
    --pred pipeline_swissprot=pipeline/pipeline_swissprot.gff3 \
    --pred pipeline_no_ecoli=pipeline/pipeline_swissprot_no_ecoli.gff3 \
    --outdir results

step "5/5  sanity checks (self-score + 10-gene deletion)"
python3 sanity_checks.py

echo; echo "Done. Baseline table:"; cat "$P/results/metrics.tsv"
