# What the numbers in the meeting deck mean

A slide-by-slide guide to every number in `meeting_ppt_v2.pptx`: what it measures, how it was calculated, and how to explain it in one sentence.

All numbers come from the final validation run (after the naming fix, the clan overlap rule and the combiner fix). The source files are in `validation/results/ecoli/`.

---

## Key terms

| Term | Meaning | Formula |
|---|---|---|
| **Sensitivity** (recall) | Of the real things, how many did we find? | found ÷ all real |
| **Precision** | Of the things we predicted, how many are real? | correct ÷ all predicted |
| **True positive (TP)** | We predicted it and it's in the ground truth | |
| **False positive (FP)** | We predicted it but it's not in the ground truth | |
| **False negative (FN)** | It's in the ground truth but we missed it | |
| **Ground truth** | The official reference we compare against | RefSeq (genes, names), UniProt (Pfam), GO Consortium (GO) |
| **Full Swiss-Prot** | Homology search against all 575,748 Swiss-Prot proteins, including E. coli's own | the "known genome" case |
| **Without E. coli relatives** | Same search, with all 30,225 *Escherichia* and *Shigella* proteins removed (545,523 left) | simulates a new, unstudied organism |

---

## Slides 1–4: the pipeline

| Number | Meaning |
|---|---|
| **4.64 Mb** (4,641,652 bp) | Size of the E. coli K-12 MG1655 genome (one chromosome) |
| **4,319 genes + proteins** | Genes Pyrodigal predicted. Each one is also translated into a protein |
| **6,601 domains on 4,119 proteins** | Pfam domains found after the clan overlap rule; 4,119 of the 4,319 proteins have at least one domain |
| **GO for 2,366 proteins** | Proteins that received at least one GO term through pfam2go |
| **40% identity, 80% coverage** (slide 4) | The rule for copying a Swiss-Prot name: at least 40% identical amino acids, and the alignment covers at least 80% of both proteins |

**How to explain the four homology statuses (slide 4):**

| Status | What happened |
|---|---|
| `confident` | A good match with a real functional name, so the name is copied |
| `uncharacterized_match` | A good match, but every good match is named "Uncharacterized…", "UPF0…" or just "Protein Yxxx" |
| `weak_match` | Something similar, but below the 40/80/80 thresholds, so it's called "hypothetical protein" |
| `no_match` | DIAMOND found nothing, so it's called "hypothetical protein" |

---

## Slide 6: validation design

Every predicted gene is first matched to a RefSeq gene **by its stop codon** (same contig, strand and stop position). The stop is fixed; the start is the uncertain part. The matched gene's UniProt ID then links it to the Pfam and GO ground truth.

---

## Slide 7: step 1, gene prediction

| Number | Meaning | How it was calculated |
|---|---|---|
| **4,287** | Real genes used as ground truth | RefSeq has 4,318 protein-coding records. Excluded: 18 pseudogenes and 3 genes that need a frameshift to read (`prfB`, `dnaX`, `copA`), leaving 4,297. A further 19 share a stop codon with another record (same gene with alternative starts, or overlapping genes), and are counted as 9, giving **4,287** |
| **3,889 exact** | Predicted with the correct start and stop | |
| **229 wrong start** | Correct stop, but a different start | |
| **169 missed** | Real genes with no prediction | |
| **201 extra** | Predictions with no real gene | |
| **96.06% sensitivity** | Share of real genes found | (3,889 + 229) ÷ 4,287 |
| **95.35% precision** | Share of predictions that are real genes | (3,889 + 229) ÷ 4,319 |
| **94.44% exact-start rate** | Share of found genes with the exact start | 3,889 ÷ (3,889 + 229) |

**One sentence:** Pyrodigal finds 96% of genes, and almost all its predictions are real; the errors are mostly small genes and start positions.

---

## Slide 8: where gene-prediction errors occur

| Number | Meaning |
|---|---|
| **24%** of genes under 150 bp found (n = 155) | The main weakness: tiny genes such as leader peptides |
| **89%** at 150–299 bp; **99–100%** above 300 bp | Detection is nearly perfect once genes are a normal size |
| **89 too long, 140 too short** | Of the 229 wrong starts, 89 predictions start too early (protein too long) and 140 start too late (protein too short) |

**One sentence:** missed genes are almost all very short, and start errors cluster close to the true start.

---

## Slide 9: a step 1 error that affects step 3

| Number | Meaning |
|---|---|
| **5** domain clans missed | Pfam domains UniProt lists but we missed, counted at clan level. All 5 are on genes predicted with a start that is too late: `fdnG`, `fdoG`, `fdhF` (selenoproteins) plus `stfR` and `insK` |
| **468–1,365 bp** | How far too late those starts are (156–455 amino acids lost), enough to remove a whole domain |
| **97.30%** | Pfam recall for the **139** proteins whose start is too late |
| **99.78%** | Pfam recall for the **3,879** proteins with exact starts |

**Why selenoproteins:** selenocysteine is encoded by TGA, normally a stop codon. Pyrodigal treats it as a stop, so it predicts a gene that starts after it.

**One sentence:** a gene-prediction error can silently remove a domain later in the pipeline.

---

## Slides 10–11: step 2, homology names

Calculated on the **4,118** genes matched to RefSeq in step 1 (3,889 + 229).

**Categories (slide 10):**

| Category | Rule |
|---|---|
| Correct: same gene | The Swiss-Prot match has the same gene symbol as RefSeq (ignoring copy numbers, e.g. `insH10` = `insH11`) |
| Correct: similar name | Different symbol, but at least half the meaningful words in the names are shared |
| Different name (review) | Both give a function, but the names don't match automatically; needs a human check |
| Function missed | RefSeq knows the function; we gave no functional name |
| Function claimed, unknown | We gave a functional name; RefSeq says the function is unknown |
| Correctly unknown | Both say unknown |

**Numbers (slide 11):**

| | Full Swiss-Prot | Without E. coli relatives |
|---|---|---|
| Correct: same gene | 3,255 | 1,838 |
| Correct: similar name | 80 | 222 |
| Different name (review) | 80 | 318 |
| Function missed | 351 | 1,388 |
| Function claimed, unknown | 21 | 9 |
| Correctly unknown | 331 | 343 |
| **Known functions recovered** | **88.56%** | **54.70%** |

**Known functions recovered** = correct ÷ genes whose function RefSeq knows:

- Full: (3,255 + 80) ÷ (3,255 + 80 + 80 + 351) = 3,335 ÷ 3,766 = **88.56%**
- Without relatives: (1,838 + 222) ÷ (1,838 + 222 + 318 + 1,388) = 2,060 ÷ 3,766 = **54.70%**

**One sentence:** with E. coli in the database, homology mostly finds the genome itself; for a new organism it recovers only about half the known functions, and the loss is mostly *missing* names, not *wrong* ones.

---

## Slides 12–13: step 3, Pfam domains

Compared on **4,106** matched genes that have a UniProt entry.

| Number | Meaning |
|---|---|
| **6,155** true positives | Pfam families that both we and UniProt assign to the same protein |
| **99.69%** precision and recall | Family level: 6,155 TP, 19 FP, 19 FN, so 6,155 ÷ 6,174 for both |
| **Clan level:** 5,832 TP, 14 FP, 5 FN; **99.76%** precision, **99.91%** recall | Related families (same Pfam clan) counted as one |
| **Before vs after the overlap rule** | Before: 8,320 domains, precision 80.1%, recall 99.9%. After: 6,601 domains, precision 99.7%, recall 99.7% |
| **19 family-level false positives** | 5 belong to the same clan as a correct family (near-duplicates); 14 are genuinely different families |

Per protein:

| Outcome | Proteins |
|---|---|
| Identical domain sets | 3,973 |
| No domains on either side | 102 |
| Only missing families | 13 |
| Only extra families | 12 |
| Both extra and missing | 6 |

**Why the chart's "Family" and "Family + overlap rule" bars look the same:** the pipeline now applies the rule itself, so both are measured after it. The real "before" is the 80.1% above.

**The overlap rule:** when two domains from the same Pfam clan overlap on a protein, keep only the higher-scoring one. Example: *thrA* had both ACT_7 and ACT_9 on the same two regions; ACT_7 was dropped.

**Caveat:** UniProt's Pfam assignments come from the same Pfam models, so this shows agreement with the reference rather than independent accuracy.

---

## Slides 14–15: step 4, GO terms

Compared against the GO Consortium file for E. coli MG1655.

| Number | Meaning |
|---|---|
| **3,363** proteins | Matched genes that have at least one experimental GO annotation |
| **60.96%** coverage | Share of those proteins that got any GO term from our pipeline |
| **63.09%** precision | Of our GO terms, the share supported by experimental GO |
| **20.07%** recall | Of the experimental GO terms, the share we predicted |
| **25.55%** / **8.22%** | Precision and recall if only exact terms count, ignoring the GO hierarchy |
| **99.45%** precision | Against all GO, including computer-inferred annotations |

By category, against experimental GO:

| Category | Precision | Recall |
|---|---|---|
| Molecular function | 64.76% | 26.70% |
| Biological process | 77.76% | 18.60% |
| Cellular component | 96.08% | 12.13% |

**How to explain:**

- **Experimental** = backed by a lab experiment (codes EXP, IDA, IPI, IMP, IGI, IEP and high-throughput versions). This is the honest benchmark.
- **The 99.45% is inflated.** Much computer-inferred GO comes from the same InterPro2GO mapping our `pfam2go` file uses, so it partly compares the method with itself.
- **The hierarchy matters.** Before comparing, every term is expanded with its parent terms. A correct but general prediction ("oxidoreductase activity" for a specific dehydrogenase) then counts as correct at the level it states. That's why precision rises from 25.55% (exact terms only) to 63.09%.
- **63% is an underestimate.** Experimental GO is incomplete: a correct prediction counts as wrong if nobody has tested that function yet.
- **"Protein binding" (GO:0005515) is excluded** from the truth, as is standard, because it is so generic.

**One sentence:** GO terms from domains are mostly right but cover only about a fifth of what is experimentally known.

---

## Slides 16–17: step 5, the final table

The final table names every gene from homology, or from its strongest Pfam domain ("X domain-containing protein") when homology gives no usable name.

| Number | Meaning |
|---|---|
| **3,766** | Matched genes whose function RefSeq knows (the denominator below) |
| **98.17%** (full) / **96.95%** (without relatives) | Known-function genes with any functional information: a correct name, a different name, or a domain-based name |
| **88.56%** / **54.70%** | Known-function genes with a fully correct name (same as step 2) |
| **63% → 97%** | Without relatives: functional information from homology only (63%) vs the final table (97%) |
| **91% → 98%** | Same comparison with full Swiss-Prot |
| **1,594** | Domain-based names among the 4,118 matched genes, without relatives. (The full output has 1,689, because unmatched genes count too) |

How 96.95% is calculated (without relatives): (2,060 correct + 318 different name + 1,273 domain names) ÷ 3,766 = 3,651 ÷ 3,766.

**One sentence:** domain names fill most of the gap homology leaves, but they are partial ("aldo/keto reductase family" rather than "methylglyoxal reductase"). That's why only 54.7% get a fully correct name.

---

## Slide 18: conclusions

A summary of the numbers above; nothing new.

---

## Slide 19: fixes made from the validation

| Fix | Numbers | Meaning |
|---|---|---|
| Naming rule | Confident names 3,630 → 3,477 (full Swiss-Prot) | 153 names like "UPF0758 protein" or "Protein YabR" no longer count as known function |
| Clan overlap rule | 8,320 → 6,601 domains; precision 80.1% → 99.7%; recall 99.9% → 99.7% | 1,719 near-duplicate domains removed, at the cost of 14 correct families |
| Combiner fallback | Hypothetical known-function genes 335 → 69 (full Swiss-Prot) | Genes whose only homology match was "Uncharacterized…" now get a Pfam domain name if they have a known domain |

---

## Slide 20: where an agent fits

| Number | Meaning |
|---|---|
| **1,273** | Genes, without relatives, that get only a partial domain-based name although RefSeq knows the full function |
| **318** | Genes whose homology name doesn't match RefSeq automatically (the "different name" cases) |
| **115** | Known-function genes still left hypothetical, without relatives |
| **yhdZ → AapP** | Example of a borrowed name: *yhdZ*'s best informative match was AapP from *Rhizobium* (71.6% identity). The top matches agree on "amino-acid ABC transporter" but disagree on the substrate, so the borrowed name is more specific than the evidence supports |

---

## Slide 21: next steps

| Number | Meaning |
|---|---|
| **Bakta v1.12.2, database v6.0 full** | The Bakta run already completed on the same genome. Its own summary reports 4,295 CDS, 78 small ORFs and **32 hypothetical proteins**. It hasn't been put through this validation yet |

---

## Likely questions, with short answers

- **Why not use BLAST?** DIAMOND does the same kind of protein search, typically hundreds of times faster with similar sensitivity. All 4,319 proteins took about 30 seconds.
- **Why remove E. coli relatives?** Swiss-Prot contains E. coli's own proteins, so the full run mostly finds the genome itself. Removing *Escherichia* and *Shigella* shows how the pipeline would do on a new organism.
- **Is 99.7% Pfam agreement impressive?** Partly. UniProt uses the same Pfam models, so it mainly confirms the search is set up correctly.
- **Why is GO recall so low?** pfam2go only links domains to broad terms, and many domains have no GO link at all.
- **What is still unvalidated?** The "different name" cases and the domain-based names need the manual reviews (`step2_manual_review.tsv`, `step5_domain_name_review.tsv`).
