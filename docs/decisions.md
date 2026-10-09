# Design decisions log

Format: date, decision, alternatives considered, reason. These become the
"justification of design decisions" in the dissertation.

## 2026-10-09: Two agents, CAM and PAM (approved by Kevin Bryson)
- **Decision:** CAM does gene detection and translation; PAM does functional annotation, mirroring AGMIAL (Bryson et al., 2006).
- **Alternatives:** Annotator + Curator agents; a single agent.
- **Reason:** Follows a proven annotation strategy and gives each agent a clear, separately evaluable job.

## 2026-10-09: Pyrodigal for gene calling
- **Alternatives:** Prodigal binary, GeneMarkS-2, Glimmer.
- **Reason:** Same algorithm as Prodigal (used by Prokka and Bakta) but a Python library, so it runs inside a graph node without subprocess calls.

## 2026-10-09: Translation cross-checked with Biopython
- **Reason:** Two independent translations of each complete CDS; any disagreement is flagged, not silently fixed.

## 2026-10-09: CAM v0 has no LLM
- **Reason:** Gene calling is well solved by deterministic tools. The LLM reasoning belongs where evidence conflicts: in PAM and in the CAM ⇄ PAM loop.

## 2026-10-09: Self-training threshold of 100 kb
- **Reason:** Pyrodigal warns below 100 kb; shorter inputs use metagenomic mode.

## 2026-10-09: Evaluation matches genes by stop codon
- **Reason:** Standard in gene-finder benchmarks (Prodigal paper). It separates "found the gene" from "found the right start".

## 2026-10-09: Test genomes are L. sakei 23K and E. coli K-12
- **Reason:** L. sakei was annotated with AGMIAL; E. coli K-12 has the best-curated reference annotation (EcoCyc).
