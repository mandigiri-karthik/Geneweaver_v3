# Phase 0 baseline: deterministic pipeline on E. coli K-12 MG1655

GFF3 exports of the deterministic pipeline (Pyrodigal -> DIAMOND/Swiss-Prot -> PyHMMER/Pfam -> pfam2go
-> combiner) for Phase 0 scoring. Do not edit these files.

| File | Homology database | Purpose |
|---|---|---|
| `pipeline_swissprot.gff3` | full Swiss-Prot (575,748 proteins) | "known genome" case, leaky like Bakta |
| `pipeline_swissprot_no_ecoli.gff3` | Swiss-Prot without *Escherichia*/*Shigella* (545,523) | fair "novel organism" case |

Both use the same gene predictions, so CDS coordinates are identical; only the functional fields differ.

## Pipeline run

| Field | Value |
|---|---|
| Run | 2026-10-16 01:10 (local), outputs in `Structural_annotation/*_results/` |
| Code | `main.py` + `Structural_annotation/` on branch `dev`, committed as `4d6fd7d` (2026-10-16 11:29) |
| Genome | NC_000913.3, md5 `92c997bcd88e983ffdb21b2712ed3736` (same bytes as `../reference/`) |

The run happened before the code was committed, so `4d6fd7d` is the closest recorded code version.

## Tool versions

| Tool | Version |
|---|---|
| Pyrodigal | 3.7.1 (single mode, trained on the genome, translation table 11) |
| DIAMOND | 2.2.8 (`blastp --more-sensitive --evalue 1e-5 --max-target-seqs 5`) |
| PyHMMER | 0.12.3 (Pfam-A.hmm is HMMER3/f format) |
| Biopython | 1.88 |

Naming rule: copy a Swiss-Prot name when identity >= 40% and both query and subject coverage >= 80%.

## Databases

Release numbers are not stored in these files, so each is identified by md5 and file date.

| File | Date | md5 |
|---|---|---|
| `uniprot_sprot.fasta.gz` | 2026-10-14 | `bc9d398533e6df582b563c6c03093bd0` |
| `uniprot_sprot_no_ecoli.fasta.gz` | 2026-10-15 | `cf527fd5843f613d9915f0d194f7ade9` |
| `swissprot.dmnd` | 2026-10-14 | `f4174be627758fbe929cb9f514a6d5f3` |
| `swissprot_no_ecoli.dmnd` | 2026-10-15 | `54b3948efcc98ee0eb3a71718e7db2e2` |
| `Pfam-A.hmm` | 2026-01-22 | `589ccc8370eaf079022db0dbb25c8aa2` |
| `Pfam-A.clans.tsv.gz` | 2026-10-16 | `a70eeffcd8105ff094ea5834456d8bab` |
| `pfam2go` (version date 2026/07/06) | 2026-10-14 | `7d9f1c4e1f3f4b811a7fec9c4d7f2c27` |

All in `Structural_annotation/databases/` (gitignored; too large for Git).

## Export

Run from `Phase-0/`:

```bash
S=../Structural_annotation
G=reference/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna

python3 export_pipeline_gff3.py --annotation pipeline/source/ecoli_swissprot_annotation.tsv \
    --genes-gff pipeline/source/ecoli_genes.gff --swissprot $S/databases/uniprot_sprot.fasta.gz \
    --genome $G --output pipeline/pipeline_swissprot.gff3

python3 export_pipeline_gff3.py --annotation pipeline/source/ecoli_swissprot_no_ecoli_annotation.tsv \
    --genes-gff pipeline/source/ecoli_genes.gff --swissprot $S/databases/uniprot_sprot_no_ecoli.fasta.gz \
    --genome $G --output pipeline/pipeline_swissprot_no_ecoli.gff3
```

`source/` holds frozen copies of the pipeline's output tables (from `Structural_annotation/`), so
the export can be repeated even after the pipeline is re-run.

### What goes into each CDS line

| GFF3 field | Taken from |
|---|---|
| seqid, start, end, strand | combined table (= Pyrodigal GFF) |
| score | Pyrodigal score |
| `product` | the pipeline's chosen product (homology, Pfam domain name, or "hypothetical protein") |
| `gene` | `GN=` of the matched Swiss-Prot entry, **only when `product_source=homology`**; Pfam-named and hypothetical genes get no symbol |
| `Dbxref` | `Pfam:` accessions and `GO:` terms |
| `product_source`, `homology_status`, `homology_match`, `homology_identity` | combined table |
| `partial`, `start_type` | Pyrodigal GFF |

**No EC numbers.** The pipeline does not assign EC numbers (Swiss-Prot FASTA headers do not carry them),
so the pipeline scores 0 on EC agreement by construction.

## Counts

| | full Swiss-Prot | no E. coli |
|---|---|---|
| CDS | 4,319 | 4,319 |
| product from homology | 3,767 | 2,444 |
| product from Pfam domain | 466 | 1,689 |
| hypothetical (no source) | 86 | 186 |
| with gene symbol | 3,767 | 2,427 (17 matched entries have no GN) |
| with Pfam/GO Dbxref | 4,119 | 4,119 |

## Validation

Both GFF3 files were read back: every line has 9 columns, every product (after percent-decoding) and
coordinate matches the source table, and coordinates match the Pyrodigal GFF line by line.

## Checksums

| File | md5 |
|---|---|
| `pipeline_swissprot.gff3` | `282aff6cb9bec3a447de8f84eeb98890` |
| `pipeline_swissprot_no_ecoli.gff3` | `b0ea0ed4736a3fa05c69dbd9c365f432` |
| `source/ecoli_swissprot_annotation.tsv` | `8478213c4ee7a5de9604f65a3c841af9` |
| `source/ecoli_swissprot_no_ecoli_annotation.tsv` | `58f3e8f9100ab20b04d82bb61ef2c58a` |
| `source/ecoli_genes.gff` | `124dd8a3be200b4184869c46df3fd9e0` |
