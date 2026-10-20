# Phase 0 reproducible bundle

Runs the Phase 0 analysis — deterministic pipeline, GFF3 export, scoring against the frozen
reference, and the scorer's own sanity checks — identically on a laptop, the SOCS OpenShift
cluster, or any machine with a container runtime. Databases are fetched and md5-verified, not
copied around.

Drop this directory at the repository root (`Geneweaver_v3/phase0-bundle/`).

## Before anything else: two repository fixes

1. **Merge `dev` into `main`.** `homology_annotator.py`, `pfam_annotator.py`, `go_mapper.py` and
   `annotation_combiner.py` live only on `dev` (`a66b833`, `4d6fd7d`). `run_phase0.sh` checks for
   them and stops with a clear message if they are missing.
2. **Apply `patches/0001-diamond-on-path.patch`.** `HomologyAnnotator` currently shells out via
   `conda run -n diamond`, which does not exist in a container. The patch makes the conda wrapper
   optional (`conda_env=None` by default) and calls `diamond` on `PATH`:

   ```bash
   git apply phase0-bundle/patches/0001-diamond-on-path.patch
   ```

   Behaviour is unchanged when you do pass `--conda-env diamond`.

Then commit `Phase-0/` and `validation/` — both are currently untracked, and they hold the frozen
reference, the scorer and every metric.

## Quick start

```bash
cd phase0-bundle
make image                       # build the container (~5 min)
make fetch   DB=/path/to/databases   # download + verify ~2.8 GB, once
make run     DB=/path/to/databases   # full Phase 0
```

Without a container (needs `diamond` on `PATH`):

```bash
make local-fetch DB=/path/to/databases
make local-run   DB=/path/to/databases
```

## Running in JupyterLab on the IDA / SOCS OpenShift cluster

Use `cluster/` rather than the container. A Jupyter pod has no root, so nothing here needs it:
the Python packages go in with `pip --user` and DIAMOND is a static binary dropped into a user
prefix.

```
cluster/
├── bootstrap.sh          no-root setup: pinned pip packages + DIAMOND 2.2.8 into $PREFIX/bin
├── paths.example.json    copy to paths.json; where the repo, databases and Bakta DB live
├── check_env.py          verifies tools, repo files and all 7 databases; exits non-zero on problems
├── run_bakta_local.sh    reproduce the frozen Bakta comparator with a local install
├── 00_setup.ipynb        drives bootstrap.sh + check_env.py
└── 01_phase0.ipynb       the five-stage run, with the metrics table at the end
```

Order: run `00_setup.ipynb`, edit the `cluster/paths.json` it creates, re-run the verify cell until
it passes, then `01_phase0.ipynb`.

If `$HOME` is small, install onto a PVC: `PREFIX=/data/tools bash cluster/bootstrap.sh`.

**The notebooks hold no logic.** Every stage shells out to the same script a terminal run would
use, so the notebook and the command line cannot drift apart, and the notebook is disposable.
The final cell compares `function agrees` against the committed baseline (90.74 / 87.71 / 51.65 %)
and prints `CHANGED` if any value moved — treat that as a finding to explain, not a number to
adopt.

### Regenerating the Bakta comparator locally

`Phase-0/bakta/` came from the **web API** (backend 1.1.0), pinned to Bakta v1.12.2 / DB v6.0 full.
`cluster/run_bakta_local.sh` runs a local install with flags matching that job config and writes to
`Phase-0/bakta_local/` — never into the frozen directory. Two flags are load-bearing:
`--keep-contig-headers` (without it Bakta renames contigs and `score_gff.py`, which matches on
sequence id, reports everything as missed) and `--replicons` (declares NC_000913.3 circular, which
affects features at the origin). The script prints the command to score both Bakta runs side by
side; do not swap the comparator until they agree.

## On OpenShift as a batch Job

The image is built root-last and runs as a non-root UID with group-writable `/work` and `/data`,
which is what OpenShift's arbitrary-UID security context requires. Push the image to the internal
registry, mount a PersistentVolumeClaim at `/data` for the databases and the repository at `/work`,
then run `scripts/fetch_databases.sh` once as a Job and `scripts/run_phase0.sh` as the analysis Job.

Note the scoring step itself needs no cluster at all: `score_gff.py` is pure standard library and
grades 4,300 genes in seconds. The cluster buys you a pinned environment and a place to keep the
databases — the GPU nodes matter for arm C, not for this.

## Contents

| Path | Purpose |
|---|---|
| `Containerfile` | Python 3.12 + DIAMOND 2.2.8 + pinned pip packages |
| `requirements.lock.txt` | exact versions from `Phase-0/pipeline/PROVENANCE.md` |
| `manifest.tsv` | database md5s, all seven verified against the originals on 2026-10-20 |
| `scripts/fetch_databases.sh` | download, verify, derive the no-*E. coli* copy, build DIAMOND indexes |
| `scripts/make_no_ecoli.py` | removes *Escherichia*/*Shigella* from Swiss-Prot (stdlib only) |
| `scripts/run_pipeline.py` | **the missing driver** — instantiates the four annotator classes |
| `scripts/run_phase0.sh` | the five-step end-to-end run |
| `patches/` | the DIAMOND-on-PATH fix |
| `Makefile` | the commands above |

## Why `run_pipeline.py` exists

The four annotator modules are library classes with no `__main__` block. Whatever drove them was
`Structural_annotation/stage0.ipynb`, deleted in `d0789ef` — so there was no way to re-run arm B
from a clean checkout. The driver instantiates them in order and reproduces exactly the filenames
`export_pipeline_gff3.py` expects:

```
ecoli_<db>_hits.tsv   ecoli_<db>_products.tsv   ecoli_pfam_hits.tsv
ecoli_go_terms.tsv    ecoli_<db>_annotation.tsv
```

Domain and GO evidence are computed once (they do not depend on the homology database); homology
and combination run once per leakage condition. Defaults match the recorded run: identity ≥ 40 %,
both coverages ≥ 80 %, E-value 1e-5, `--more-sensitive`, `--max-target-seqs 5`.

## Database pinning

`fetch_databases.sh` verifies every download against `manifest.tsv` and **refuses to continue on a
mismatch**, because every Phase 0 number is tied to these exact files.

- **Swiss-Prot** is pinned to UniProt release **2026_03**. That is still `current_release` as of
  2026-10-20, so the plain URL works today. When 2026_04 ships, the md5 check will fail and the
  script prints the `previous_releases/release-2026_03/` URL to use instead. This is the most
  likely thing to break; it fails loudly rather than silently changing your baseline.
- **The no-*E. coli* copy is derived, not downloaded.** `make_no_ecoli.py` was validated against
  the original: 575,748 entries in, **545,523 kept, 30,225 removed, and the kept accession set is
  identical to `uniprot_sprot_no_ecoli.fasta.gz`**. It asserts the 545,523 count and fails if it
  differs. Note gzip byte-identity depends on the zlib build, so the md5 may differ even when the
  entry set is exact — the counts are the real check.
- **Pfam** release is a variable (`PFAM_RELEASE`, default `Pfam38.2`). Your `Pfam-A.hmm` is dated
  2026-01-22 with md5 `589ccc8370eaf079022db0dbb25c8aa2`; I could not determine which release that
  is without downloading several multi-GB candidates. If the md5 check fails, set `PFAM_RELEASE` to
  another release from <https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/> and re-run.
- **`.dmnd` indexes are rebuilt locally**, never downloaded — their bytes depend on the DIAMOND
  build. The recorded md5s are informational and are not verified.
- **pfam2go** comes from the EBI GO mirror (`ftp.ebi.ac.uk/pub/databases/GO/goa/external2go/pfam2go`)
  rather than `current.geneontology.org`.

## What was and was not tested

Verified here:

- all seven recorded database md5s match the files in `Structural_annotation/databases/`
- `make_no_ecoli.py` reproduces the leakage-control entry set **exactly**
- every download URL returns HTTP 200, and DIAMOND v2.2.8 exists as a Linux release
- UniProt `current_release` is still 2026_03
- all scripts pass syntax checks; both Python entrypoints have working `--help`

**Not tested:** the container build itself — no container runtime was available in the environment
where this was written. Run `make image` once and expect to fix small things (most likely the
`diamond` tarball layout). The pipeline has also not been executed end to end here, because the
`dev`-branch modules and the 2.8 GB of databases were not both present in one place.
