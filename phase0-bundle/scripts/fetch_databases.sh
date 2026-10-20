#!/usr/bin/env bash
# Download and verify the Phase 0 reference databases. Idempotent: skips files
# that already verify. Writes into $GENEWEAVER_DB (default /data).
set -euo pipefail

DB="${GENEWEAVER_DB:-/data}"
UNIPROT_RELEASE="${UNIPROT_RELEASE:-2026_03}"
PFAM_RELEASE="${PFAM_RELEASE:-Pfam38.2}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="$HERE/../manifest.tsv"

mkdir -p "$DB"; cd "$DB"

md5_of() { if command -v md5sum >/dev/null; then md5sum "$1" | cut -d" " -f1; else md5 -q "$1"; fi; }

expected() { awk -F"\t" -v f="$1" '$1==f{print $2}' "$MANIFEST"; }

verify() {  # verify FILE -> 0 if md5 matches manifest
    local f="$1" want got
    want="$(expected "$f")"; [ -z "$want" ] && return 1
    [ -f "$f" ] || return 1
    got="$(md5_of "$f")"
    [ "$got" = "$want" ] && return 0
    echo "  md5 MISMATCH for $f" >&2
    echo "    expected $want" >&2
    echo "    got      $got" >&2
    return 1
}

fetch() {  # fetch URL FILE
    local url="$1" f="$2"
    if verify "$f" 2>/dev/null; then echo "ok (cached)  $f"; return; fi
    echo "downloading  $f"
    curl -fSL --no-progress-meter --retry 3 --retry-delay 5 -o "$f.part" "$url"
    mv "$f.part" "$f"
}

echo "== Swiss-Prot (UniProt $UNIPROT_RELEASE) =="
# current_release serves whatever is newest; once 2026_04 ships, the md5 below
# will stop matching and you must switch to the previous_releases path.
SP_URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
fetch "$SP_URL" uniprot_sprot.fasta.gz
if ! verify uniprot_sprot.fasta.gz; then
    cat >&2 <<EOF

  The Swiss-Prot file does not match the md5 recorded for release $UNIPROT_RELEASE.
  UniProt has almost certainly published a newer release. Re-fetch the pinned one:

    curl -fSL -o $DB/uniprot_sprot.fasta.gz \\
      https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-$UNIPROT_RELEASE/knowledgebase/complete/uniprot_sprot.fasta.gz

  Do NOT proceed with a different release: every number in Phase 0 is tied to this one.
EOF
    exit 1
fi

echo "== leakage-control copy (Escherichia / Shigella removed) =="
# The md5 cannot be relied on here: gzip output varies with the zlib build, so the
# derived file legitimately differs byte-for-byte from the Phase 0 original while
# holding exactly the same entries. The entry count is the real integrity check.
count_entries() { python3 -c "
import gzip,sys
print(sum(1 for l in gzip.open(sys.argv[1],'rt') if l.startswith('>')))" "$1" 2>/dev/null || echo 0; }

if [ -f uniprot_sprot_no_ecoli.fasta.gz ] && [ "$(count_entries uniprot_sprot_no_ecoli.fasta.gz)" = "545523" ]; then
    echo "ok (cached)  uniprot_sprot_no_ecoli.fasta.gz  (545523 entries)"
else
    python3 "$HERE/make_no_ecoli.py" \
        --input uniprot_sprot.fasta.gz \
        --output uniprot_sprot_no_ecoli.fasta.gz \
        --expect-kept 545523
    # make_no_ecoli.py already asserts --expect-kept 545523 and exits non-zero otherwise.
    echo "   derived and entry count asserted (545523)"
fi

echo "== Pfam ($PFAM_RELEASE) =="
if verify Pfam-A.hmm 2>/dev/null; then
    echo "ok (cached)  Pfam-A.hmm"
else
    fetch "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/$PFAM_RELEASE/Pfam-A.hmm.gz" Pfam-A.hmm.gz
    gunzip -kf Pfam-A.hmm.gz
    verify Pfam-A.hmm || {
        echo "  Pfam-A.hmm does not match the recorded md5." >&2
        echo "  Phase 0 used a file dated 2026-01-22; set PFAM_RELEASE to the matching" >&2
        echo "  release and re-run. Available: https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/" >&2
        exit 1
    }
fi
fetch "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/$PFAM_RELEASE/Pfam-A.clans.tsv.gz" Pfam-A.clans.tsv.gz

echo "== pfam2go =="
# pfam2go is a ROLLING file: the EBI/GO mirror always serves the newest mapping and
# keeps no versioned archive, so downloading it cannot reproduce the baseline. The
# exact file the Phase 0 numbers were produced with (GO version date 2026/07/06) is
# vendored in the repo at Phase-0/reference/pfam2go - only 649 KB - so copy that.
VENDORED="$HERE/../../Phase-0/reference/pfam2go"
if verify pfam2go 2>/dev/null; then
    echo "ok (cached)  pfam2go"
elif [ -f "$VENDORED" ]; then
    cp "$VENDORED" pfam2go
    verify pfam2go && echo "ok           pfam2go (from Phase-0/reference/)" || {
        echo "  the vendored pfam2go does not match the manifest md5" >&2; exit 1; }
else
    echo "  Phase-0/reference/pfam2go is missing; falling back to the live mirror." >&2
    echo "  WARNING: the live file is newer than the baseline and will change GO results." >&2
    fetch "https://ftp.ebi.ac.uk/pub/databases/GO/goa/external2go/pfam2go" pfam2go
    verify pfam2go || echo "  (expected - the live mapping has moved on since 2026/07/06)" >&2
fi

echo "== DIAMOND indexes (rebuilt locally; md5 depends on the DIAMOND build) =="
for base in swissprot swissprot_no_ecoli; do
    src="uniprot_sprot.fasta.gz"
    [ "$base" = "swissprot_no_ecoli" ] && src="uniprot_sprot_no_ecoli.fasta.gz"
    if [ -f "$base.dmnd" ] && [ ! "$src" -nt "$base.dmnd" ]; then
        echo "ok (cached)  $base.dmnd"
    else
        [ -f "$base.dmnd" ] && echo "   $src is newer than $base.dmnd - rebuilding"
        echo "building     $base.dmnd"
        diamond makedb --in "$src" --db "$base" --quiet
    fi
done

echo
echo "All databases present in $DB"
ls -lh "$DB"
