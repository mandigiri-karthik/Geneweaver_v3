# Hand check of the scoring rules

Goal: confirm, by eye and without an LLM, that the scorer's fixed string rules give the right answer.
Each tool has a sheet of **50 random matched genes** (seed 42), produced by `score_gff.py`:

| Sheet | Suggested order |
|---|---|
| `pipeline_no_ecoli_handcheck_sample.tsv` | 1st: the fair baseline that later phases must beat; 18/50 are `wrong_function` |
| `bakta_handcheck_sample.tsv` | 2nd |
| `pipeline_swissprot_handcheck_sample.tsv` | 3rd, optional (same genes as no-E. coli, mostly the same calls) |

Open a sheet in Excel or Numbers (tab-separated). For background on a gene, look up its `locus_tag`
(b-number) in EcoCyc (ecocyc.org) or UniProt.

## Columns

| Column | Meaning |
|---|---|
| `locus_tag`, `strand`, `stop` | the reference gene (stop is the shared coordinate that made it a match) |
| `ref_start` / `pred_start` | start coordinates; `start_exact` says whether they are equal |
| `ref_symbols` / `pred_gene` | reference gene symbols (RefSeq + UniProt, lowercased) / predicted symbol |
| `ref_ec` / `pred_ec` | reference EC numbers (UniProt) / predicted ECs |
| `ref_product` / `pred_product` | RefSeq product / predicted product |
| `symbol_match`, `ec_match`, `product_match` | the scorer's rule results (`n/a` = no full reference EC) |
| `category` | the scorer's verdict for this gene |

## Fill in (Y or N)

| Column | Question |
|---|---|
| `hand_same_gene` | Is the predicted CDS the same gene as the reference (same locus, sensible start)? |
| `hand_same_function` | Does the predicted name describe the same function as the reference, even if worded differently? A vague domain name that is compatible but not specific (e.g. "TrkA-N domain-containing protein" for K⁺ channel Kch) is **N** — write "vague" in the note. |
| `hand_category_ok` | Is the scorer's `category` right? |
| `hand_note` | anything notable: "vague", "paralog", "naming convention", "reference outdated"… |

## Summarise

```bash
cd Phase-0
python3 summarize_handcheck.py results/*_handcheck_sample.tsv
```

It prints, per sheet, the share of Y answers for each hand column, and lists the genes the scorer
called `wrong_function` where you said the function is the same. That list is the evidence of how
strict the string rules are, and the share of `hand_category_ok = Y` is the number to report.

## What is already known

The rules are deliberately strict, so the agreement numbers are a **lower bound**:
- Product names must match after normalisation; "Bifunctional aspartokinase/homoserine dehydrogenase 1"
  does not match "fused aspartate kinase/homoserine dehydrogenase 1".
- Pfam-domain names from the pipeline ("EAL domain-containing protein") usually share no symbol, EC or
  product with the reference ("putative cyclic di-GMP phosphodiesterase PdeD") and land in `wrong_function`.
- IS element proteins are named differently by Bakta ("IS5 family IS5 transposase", `tnp`) and RefSeq
  ("IS5 transposase and trans-activator", `insH21`).
