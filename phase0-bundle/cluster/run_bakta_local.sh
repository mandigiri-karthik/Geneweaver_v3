#!/usr/bin/env bash
# Reproduce the frozen Phase 0 Bakta comparator with the local Bakta install on the IDA pod.
#
# The frozen result in Phase-0/bakta/ came from the web API (backend 1.1.0, Bakta v1.12.2,
# DB v6.0 full) with job config {completeGenome, translationTable 11, keepContigHeaders,
# minContigLength 1, compliant false}. The flags below match the command the server ran.
# Two of them are load-bearing:
#
#   --keep-contig-headers  without it Bakta renames contigs, and score_gff.py matches on
#                          sequence id, so scoring reports everything as missed for a reason
#                          that has nothing to do with annotation quality
#   --replicons            declares NC_000913.3 circular; affects features at the origin,
#                          including the truncated CDS at 3..98
#
# Output goes to Phase-0/bakta_local/ - never into Phase-0/bakta/, which is frozen.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PATHS="${PATHS:-$HERE/paths.json}"

eval "$(python3 - "$PATHS" <<'EOF'
import json, sys
c = json.load(open(sys.argv[1]))
for k, d in (("repo", ""), ("bakta_db", ""), ("bakta_env", ""), ("threads", 4)):
    print(f'{k.upper()}="{c.get(k, d)}"')
EOF
)"

[ -n "$BAKTA_DB" ]  || { echo "ERROR: set \"bakta_db\" in $PATHS"  >&2; exit 2; }
[ -n "$BAKTA_ENV" ] || { echo "ERROR: set \"bakta_env\" in $PATHS" >&2; exit 2; }

OUTDIR="${OUTDIR:-$REPO/Phase-0/bakta_local}"
GENOME="$REPO/Phase-0/reference/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"
REPLICONS="$REPO/Phase-0/bakta/replicons.tsv"

for f in "$GENOME" "$REPLICONS"; do
    [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 2; }
done
[ -d "$BAKTA_DB" ] || { echo "ERROR: bakta db not a directory: $BAKTA_DB" >&2; exit 2; }

# The Bakta env has no name, only a path - 'conda activate bakta' fails; use the prefix.
BAKTA=(conda run --no-capture-output -p "$BAKTA_ENV" bakta)

mkdir -p "$OUTDIR"
echo "bakta  : $("${BAKTA[@]}" --version)"
echo "env    : $BAKTA_ENV"
echo "db     : $BAKTA_DB"
echo "genome : $GENOME"
echo "out    : $OUTDIR"
echo "threads: $THREADS"
echo

"${BAKTA[@]}" \
    --db "$BAKTA_DB" \
    --prefix result \
    --output "$OUTDIR" \
    --replicons "$REPLICONS" \
    --complete \
    --keep-contig-headers \
    --translation-table 11 \
    --min-contig-length 1 \
    --gram '?' \
    --threads "$THREADS" \
    --tmp-dir "$OUTDIR/tmp" \
    --force \
    "$GENOME" 2>&1 | tee "$OUTDIR/run.log"

echo
echo "== version block from result.json (put this in PROVENANCE) =="
python3 - "$OUTDIR/result.json" <<'EOF'
import json, sys
print(json.dumps(json.load(open(sys.argv[1])).get("version", {}), indent=2))
EOF

echo
echo "== IPS count (the leakage number) =="
grep -iE 'detected IPSs?' "$OUTDIR/run.log" || echo "   not found in log"

cat <<EOF

Now compare against the frozen web-API run. Do NOT replace it until they agree:

    cd $REPO/Phase-0
    python3 score_gff.py \\
        --pred bakta_web=bakta/result.gff3 \\
        --pred bakta_local=bakta_local/result.gff3 \\
        --pred pipeline_swissprot=pipeline/pipeline_swissprot.gff3 \\
        --pred pipeline_no_ecoli=pipeline/pipeline_swissprot_no_ecoli.gff3 \\
        --outdir results_bakta_compare
EOF
