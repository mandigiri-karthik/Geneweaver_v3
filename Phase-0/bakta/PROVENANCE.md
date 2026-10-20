# Phase 0 baseline: Bakta on E. coli K-12 MG1655

Frozen Bakta output for the Phase 0 baseline. Do not edit or re-run into this folder.

## Versions

| Field | Value |
|---|---|
| Bakta | v1.12.2 |
| Database | v6.0, full |
| Web backend | 1.1.0 (api.bakta.computational.bio) |
| Citation | Schwengers et al. 2021, doi:10.1099/mgen.0.000685 |

Versions were taken from three places that agree: the API `/version` endpoint (logged in `run.log`),
the `version` block of `result.json`, and the header of `result.gff3`.

## Input

`../reference/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna`
(md5 `92c997bcd88e983ffdb21b2712ed3736`, the frozen Phase 0 reference; see `../reference/PROVENANCE.md`).

Replicon table (`replicons.tsv`): NC_000913.3, chromosome, circular.

## Exact commands

Run from `Phase-0/`:

```bash
python3 -u run_bakta.py reference/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna bakta 2>&1 | tee bakta/run.log
```

`run_bakta.py` sends this job config to the web API:

```json
{"completeGenome": true, "translationTable": 11, "keepContigHeaders": true, "minContigLength": 1, "compliant": false}
```

Bakta command run on the server (from `run.log`):

```
bakta --tmp-dir /cache --threads 8 --prefix result -o /results --db /db/db --replicons /data/replicons.tsv --complete --keep-contig-headers --gram ? /data/fastadata.fasta --force
```

## Run

| Field | Value |
|---|---|
| Job ID | 2e50d865-5033-420e-b7d6-a3e7102e6e73 |
| Submitted | 2026-10-19T16:45:33Z |
| Bakta run (server) | 2026-10-19 16:45:52 – 16:52:32 UTC (6:40) |
| Status | SUCCESSFUL |

An earlier attempt the same minute (job 3e5847d7-…) was abandoned because the client crashed on an
empty status reply; its log is kept in `run_failed_attempt1.log`. The client was then fixed to treat
that reply as pending.

The GFF3 feature coordinates (columns 1–8) are identical to the earlier Bakta run from 2026-10-15
in `Structural_annotation/bakta_results/`, so the result is reproducible with these versions.

## Annotation summary (result.txt)

| Feature | Count |
|---|---|
| CDS | 4,295 (+ 78 sORFs = 4,373 CDS rows in the GFF3) |
| hypothetical CDS | 32 |
| pseudogenes | 5 |
| tRNA / tmRNA / rRNA / ncRNA | 88 / 1 / 22 / 232 |
| ncRNA regions / CRISPR arrays / oriC | 57 / 2 / 2 |

## Notes for scoring

- **Pseudogenes are tagged `pseudogene=unitary` in the GFF3**, not `pseudo=True`. These are the 5 counted
  above (and one tRNA with `pseudogene=unknown`).
- One extra CDS (`KOIAIB_00001`, 3..98 +) carries `pseudo=True`; it is a truncated hypothetical at the
  sequence origin.
- EC numbers are in `Dbxref=EC:...`; gene symbols in `gene=`; product in `product=`.

## Leakage caveat

Bakta's database (UniRef/RefSeq-derived PSC/IPS) contains *E. coli* K-12 proteins, so many CDS are
annotated by exact identity lookup (`detected IPSs: 4246` in `run.log`). Its score on this genome is
therefore optimistic and is not a fair estimate for a novel organism.

## Files and checksums

| File | md5 |
|---|---|
| `result.gff3` | `d96302ba4a0b49fae7584f4376edc593` |
| `result.tsv` | `29280a73f382a431178b3e52c84c114b` |
| `result.faa` | `69e5658f5ff766c1904b42ce496af306` |
| `result.json` | `708d069a5e9c75594ba73c09d1cffea7` |
| `result.hypotheticals.tsv` | `ad15bbf62745896339da9a215677fd12` |
| `result.inference.tsv` | `9646faca6f3d4f9e00a6de485da5b678` |
| `result.txt` | `c30a3b1a9140f49293a4565894b8e2e0` |
