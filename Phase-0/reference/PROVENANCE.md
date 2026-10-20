# Phase 0 reference: E. coli K-12 MG1655

Frozen ground truth for every scoring run. Do not edit, replace or re-download these files.
If the reference ever has to change, make a new dated folder and re-score every tool against it.

## Source

| Field | Value |
|---|---|
| Organism | *Escherichia coli* str. K-12 substr. MG1655 (taxid 511145) |
| Assembly | GCF_000005845.2 (ASM584v2), released 2013-09-26 |
| Sequence | NC_000913.3 (= GenBank U00096.3), 4,641,652 bp, circular |
| Annotation | "Annotation submitted by NCBI RefSeq", annotation release 2026-09-02 |
| Downloaded | 2026-10-19T16:41:38Z |
| Tool | NCBI Datasets API v2 (curl) |

Download URL:

```
https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000005845.2/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT
```

Exact command:

```bash
curl -sS --fail -o ncbi_GCF_000005845.2.zip "<URL above>"
unzip -q ncbi_GCF_000005845.2.zip
```

## Files and checksums

`ncbi_GCF_000005845.2.zip` (the untouched archive): md5 `3f69280ab21096161023cb001815c4d7`

| File | md5 |
|---|---|
| `ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna` | `92c997bcd88e983ffdb21b2712ed3736` |
| `ncbi_dataset/data/GCF_000005845.2/genomic.gff` | `837628601d3bbf74c218716d7ea2df86` |
| `ncbi_dataset/data/GCF_000005845.2/protein.faa` | `4c19a285fd7a4e75ce34a466fcf0607f` |
| `ncbi_dataset/data/GCF_000005845.2/sequence_report.jsonl` | `1c42233e4ad69e6a49a6d0607bb52439` |
| `ncbi_dataset/data/assembly_data_report.jsonl` | `68e28bfd80d9509336b885ac646a3710` |
| `ncbi_dataset/data/dataset_catalog.json` | `9f4a71883c2124d969e5903e1e2dbe30` |

All six files matched NCBI's `md5sum.txt` at download. The genome FASTA is byte-identical to the
copy used in earlier runs (`E-coli_data/`), so existing Bakta and pipeline outputs used the same sequence.

Verify later with:

```bash
md5sum -c md5sum.txt      # Linux
```

## Feature counts (genomic.gff)

| Feature | Count |
|---|---|
| gene | 4,506 |
| CDS rows | 4,340 (some CDS are split over several rows) |
| pseudogene | 145 |
| ncRNA / tRNA / rRNA | 108 / 86 / 22 |
| Proteins in protein.faa | 4,300 |

NCBI's summary: 4,290 protein-coding genes, 145 pseudogenes, 215 non-coding.

## Curated, not a PGAP prediction

Evidence that this is the curated Blattner/EcoCyc annotation rather than an automatic PGAP run:

- The annotation provider is "NCBI RefSeq", not "NCBI Prokaryotic Genome Annotation Pipeline (PGAP)".
- The RefSeq record is a copy of GenBank U00096.3, which is submitted by the University of Wisconsin (Blattner lab).
- All 4,506 genes use the classic `b####` locus tags. PGAP would assign `*_RS#####` tags.
- All 4,340 CDS rows carry an `ECOCYC:` cross-reference.
- No CDS has an `inference=` qualifier, which PGAP adds to every prediction (e.g. GeneMarkS-2+).

## Known limitation

The GFF has **no EC numbers**. They come from UniProt instead (next section).

# UniProt functional reference (EC numbers, gene names)

## Source

| Field | Value |
|---|---|
| Database | UniProtKB release **2026_03** (released 2026-09-02; API deployment 2026-10-15) |
| Query | reference proteome UP000000625 (*E. coli* K-12 MG1655), all entries |
| Entries | 4,403, all reviewed (Swiss-Prot) |
| Downloaded | 2026-10-19T16:58:29Z |

Exact command (run in `reference/uniprot/`):

```bash
URL='https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000000625&format=tsv&fields=accession,reviewed,gene_primary,gene_names,gene_oln,ec,protein_name,length'
curl -sS --fail -D headers.txt -o UP000000625_ecoli_k12.tsv "$URL"
```

`headers.txt` keeps the HTTP response headers, including `x-uniprot-release`; `url.txt` holds the URL.

## Joined table: `reference_function.tsv`

Built from the frozen GFF and UniProt table (run from `Phase-0/`):

```bash
python3 build_reference_function.py \
    --gff reference/ncbi_dataset/data/GCF_000005845.2/genomic.gff \
    --uniprot reference/uniprot/UP000000625_ecoli_k12.tsv \
    --output reference/reference_function.tsv
```

One row per reference CDS (4,318 unique CDS IDs; split CDS merged into one span). Each CDS is linked
to UniProt by the `UniProtKB/Swiss-Prot` accession in its Dbxref, or, if that is missing, by its
`b####` locus tag in UniProt's ordered locus names.

| | Count |
|---|---|
| linked by Dbxref accession | 4,281 |
| linked by locus tag | 35 (e.g. alternative forms of mrcB, crl, IS element proteins) |
| unlinked | 2 (insN, insO: both pseudogene fragments) |
| CDS with >= 1 EC number | 1,742 |
| EC values in total | 1,956, of which 459 are partial (e.g. `3.2.-.-`) |
| RefSeq gene symbol differs from UniProt primary gene name | 121 (e.g. lpd/lpdA, mtn/mtnN) |
| CDS marked `pseudo=true` | 18 |

Columns: `cds_id, locus_tag, strand, start, end, pseudo, refseq_gene, refseq_product,
uniprot_accession, link_method, uniprot_gene, uniprot_protein_name, ec_numbers` (ECs `;`-separated).

Points the scorer has to decide:
- Gene symbol: the RefSeq and UniProt names disagree for 121 genes, so accept a match to either.
- Partial EC numbers (`3.2.-.-`): decide whether a predicted full EC inside a partial one counts.

## Checksums

| File | md5 |
|---|---|
| `uniprot/UP000000625_ecoli_k12.tsv` | `c73978bf3a273747ea38f4c2d4958850` |
| `reference_function.tsv` | `63c4839d5c0eaa160c7ef14e8b5815e0` |
