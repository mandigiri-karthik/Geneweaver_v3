# Phase 0: measuring tool and baseline

One scorer, `score_gff.py`, that grades any GFF3 annotation of *E. coli* K-12 MG1655 against a frozen,
curated reference. Every later phase uses this exact scorer and reference, so numbers stay comparable.

## Layout

| Path | What it is |
|---|---|
| `reference/` | frozen ground truth: NCBI RefSeq GFF3 + FASTA, UniProt 2026_03 table, `reference_function.tsv` ([PROVENANCE](reference/PROVENANCE.md)) |
| `bakta/` | frozen Bakta 1.12.2 / DB 6.0 output ([PROVENANCE](bakta/PROVENANCE.md)) |
| `pipeline/` | frozen GFF3 exports of the deterministic pipeline, two homology databases ([PROVENANCE](pipeline/PROVENANCE.md)) |
| `run_bakta.py` | runs Bakta through its web API |
| `export_pipeline_gff3.py` | turns the pipeline's combined table into GFF3 |
| `build_reference_function.py` | joins the RefSeq CDS to UniProt to add EC numbers |
| `score_gff.py` | the scorer |
| `sanity_checks.py` | self-score and deletion tests for the scorer (outputs in `sanity/`) |
| `summarize_handcheck.py` | summarises the filled-in hand-check sheets |
| `results/` | baseline table, per-gene tables, error lists, hand-check sheets |

## Running the scorer

```bash
cd Phase-0
python3 score_gff.py \
    --pred bakta=bakta/result.gff3 \
    --pred pipeline_swissprot=pipeline/pipeline_swissprot.gff3 \
    --pred pipeline_no_ecoli=pipeline/pipeline_swissprot_no_ecoli.gff3 \
    --outdir results
```

Python 3 standard library only. The reference paths default to `reference/`.

Outputs in `--outdir`:

| File | Content |
|---|---|
| `metrics.md`, `metrics.tsv` | one column/row per tool |
| `<tool>_genes.tsv` | every reference gene: matched or missed, start, function fields, category |
| `<tool>_errors.tsv` | reference genes the tool missed or called hypothetical although the function is known |
| `<tool>_false_positives.tsv` | predicted CDS that match no reference gene or pseudogene |
| `<tool>_handcheck_sample.tsv` | 50 random matched genes (seed 42), reference and prediction side by side, with empty hand-check columns |
| `errors_all_tools.tsv` | (2+ tools) every gene any tool missed or called hypothetical, all tools' outcomes side by side |

## Scoring rules

All comparisons are fixed string rules. No LLM is used as a judge.

### Reference genes

- Reference CDS come from the RefSeq GFF3. CDS split over several rows (ribosomal slippage in prfB, dnaX,
  copA) are merged into one span per CDS ID.
- A **reference gene** is one locus: sequence + strand + stop coordinate. CDS that share a stop
  (alternative starts such as mrcB, infB, clpB; 9 loci) form one gene with several valid starts.
- CDS marked `pseudo=true` are not reference genes. Their spans, together with the 145 `pseudogene`
  features, form the pseudogene regions.
- Result: **4,290 reference genes**, matching NCBI's 4,290 protein-coding genes.

### Gene recovery

- A predicted CDS **matches** a reference gene when sequence, strand and stop coordinate are the same
  (stop = end on +, start on −). Split predicted CDS are merged by ID first.
- Predicted CDS with `pseudo=True` or a `pseudogene=` attribute (Bakta's style) are not counted as
  predictions; they go into the pseudogene count.
- A prediction whose stop lies inside a reference pseudogene region on the same strand is **not** a
  false positive; it is reported separately and left out of precision.
- **Recall** = matched reference genes / reference genes.
- **Precision** = matched predicted loci / (predicted loci − predictions on reference pseudogenes).

### Start-site accuracy

Among matched genes, the share whose predicted start equals one of the reference starts.

### Function agreement (matched genes)

- **Gene symbol**: predicted `gene=` equals the RefSeq gene symbol or the UniProt primary gene name,
  ignoring case (they differ for 121 genes, e.g. lpd/lpdA). No symbol counts as no match.
- **EC number**: reference ECs come from UniProt (`reference_function.tsv`). Only **full** EC numbers
  (four numbers, e.g. `2.7.2.4`) are scored; the match is **exact string equality** of at least one
  predicted EC with one reference EC. Denominator: matched genes with at least one full reference EC.
  Genes whose reference ECs are all partial (e.g. `3.2.-.-`) are not scored and are reported as a count.
  "EC set identical" (all full ECs the same) is reported as an extra, stricter number.
  Predicted ECs are read from `Dbxref=EC:...` and `ec_number=` / `EC_number=`.
- **Product name**: both names are lowercased, every character other than letters, digits and `.`
  becomes a space, and the words `putative, probable, possible, predicted, potential, likely, protein`
  are removed. A match is equality with the normalised RefSeq product or UniProt recommended name
  (the UniProt name before the first ` (`).
- **Hypothetical**: a product is hypothetical when it is empty or every word is uninformative:
  `hypothetical, uncharacterized, putative, probable, predicted, conserved, protein, domain, containing,
  family, unknown, function, of`, a `DUF`/`UPF`/`PF` number, or a gene name like `YaaJ`. So
  "DUF2502 domain-containing protein YaaX" is hypothetical, but "Uncharacterized transporter YaaJ" is not.
  The reference **knows the function** when its RefSeq product is not hypothetical.
- Each matched gene gets one **category**:

| Category | Rule |
|---|---|
| `hypothetical_where_known` | prediction hypothetical, reference knows the function |
| `both_hypothetical` | both hypothetical |
| `named_where_reference_hypothetical` | prediction names a function, reference is hypothetical |
| `agrees` | both know a function and symbol, EC or product matches |
| `wrong_function` | both know a function and symbol, EC and product all disagree |

The rates for `agrees`, `hypothetical_where_known` and `wrong_function` use matched genes whose
reference knows the function as the denominator.

### Whole-annotation counts

- **Hypothetical share** = predicted CDS with a hypothetical product / all predicted CDS.
- **Pseudogenes** = predicted pseudogene features plus pseudo-tagged CDS (unique IDs).

## Baseline (2026-10-19)

Produced by the command under "Running the scorer"; source files in `results/`.

| Metric | Bakta | Pipeline, full Swiss-Prot | Pipeline, no *E. coli* |
|---|---:|---:|---:|
| Predicted CDS | 4,367 | 4,319 | 4,319 |
| Matched (same strand + stop) | 4,179 | 4,119 | 4,119 |
| **Recall** | **97.4%** | **96.0%** | **96.0%** |
| **Precision** | **98.7%** | **98.4%** | **98.4%** |
| False positives | 57 | 67 | 67 |
| Predictions on reference pseudogenes | 131 | 133 | 133 |
| **Start exact** (of matched) | **94.5%** | **94.4%** | **94.4%** |
| Gene symbol match (of matched) | 93.6% | 87.4% | 46.1% |
| EC exact match (of 1,357 matched genes with a full reference EC) | 97.8% | 0.0%¹ | 0.0%¹ |
| EC set identical | 91.2% | 0.0%¹ | 0.0%¹ |
| Matched genes with only partial reference EC (not scored) | 375 | 376 | 376 |
| Product match, normalised (of matched) | 47.1% | 85.4% | 44.2% |
| **Function agrees** (symbol, EC or product; of matched, ref known) | **90.7%** | **87.7%** | **51.7%** |
| Hypothetical where reference knows function | 4.8% | 2.7% | 4.3% |
| Wrong function² | 4.5% | 9.6% | 44.0% |
| Hypothetical share (of all predicted CDS) | 14.0% | 10.8% | 11.8% |
| Pseudogenes predicted | 6 | 0 | 0 |

Reference: 4,290 protein-coding genes. Function rates use matched genes whose reference function is known
(3,737 / 3,717 / 3,717).

¹ The pipeline does not assign EC numbers.
² A lower bound on quality, not a count of definite errors: the string rules cannot credit a vague but
compatible name such as "TrkA-N domain-containing protein" for the K⁺ channel Kch. Most of the
no-*E. coli* run's wrong-function genes are Pfam-domain names of this kind. The hand check
(`results/HANDCHECK.md`) measures how strict the rules are.

Bakta's product match is low (47%) because its naming style differs from RefSeq's; its symbols and ECs agree.

### Error lists (targets for Phase 3)

| File | Rows |
|---|---|
| `results/bakta_errors.tsv` | 289 (111 missed, 178 hypothetical where known) |
| `results/pipeline_swissprot_errors.tsv` | 271 (171 missed, 100 hypothetical where known) |
| `results/pipeline_no_ecoli_errors.tsv` | 332 (171 missed, 161 hypothetical where known) |
| `results/errors_all_tools.tsv` | 491 genes failed by at least one tool (146 by all three), each tool's outcome side by side |

Genes that all three tools miss are mostly short proteins (YabR, YkfO, YahV, …).

## Hand check

`results/<tool>_handcheck_sample.tsv` holds 50 random matched genes per tool, with empty Y/N columns.
Instructions are in `results/HANDCHECK.md`; summarise the filled sheets with
`python3 summarize_handcheck.py results/*_handcheck_sample.tsv`. **Status: not done yet.**

## Sanity checks

```bash
cd Phase-0
python3 sanity_checks.py          # exits 1 if any check fails; output saved in sanity/sanity_checks.log
```

1. **Self-score.** The reference GFF3 with its UniProt EC numbers added as `ec_number=` (1,746 CDS rows)
   is scored against the reference. Result: recall, precision, start, gene symbol, EC exact, EC set,
   product and function agreement are all **100%**; false positives, hypothetical-where-known and
   wrong function are **0**. The 10 duplicate stops are the reference's own alternative-start CDS.
2. **Deletion.** 10 random reference genes (seed 42: xthA, aroM, yaiS, tamA, bipA, fadL, asr, ynaL,
   casB, cueR) are removed from that copy. The scorer reports **exactly those 10** as missed
   (recall 4,280/4,290 = 99.77%) and every other metric stays at 100%.

Last run: 2026-10-19, all checks passed. The deleted genes are listed in `sanity/deleted_genes.tsv`.

## Caveat for the thesis

Bakta's database contains *E. coli* K-12 proteins (4,246 CDS were identified by exact protein identity),
and the full-Swiss-Prot pipeline run can copy names from *E. coli*'s own Swiss-Prot entries. Both scores
are therefore optimistic. `pipeline_swissprot_no_ecoli` (all *Escherichia*/*Shigella* removed) is the
fair estimate for a novel organism. This is the leakage problem Kevin raised.
