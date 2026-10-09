# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Current state

**CAM v0 is built** (gene detection + translation agent, no LLM yet):

- `src/geneweaver/cam/` — LangGraph StateGraph: `load_genome -> call_genes (Pyrodigal) ->
  translate (Biopython cross-check) -> qc_check -> [mark_for_review] -> export`.
  Typed state in `state.py`, pure node functions in `nodes.py`, graph + SqliteSaver in `graph.py`.
- `src/geneweaver/evaluation/compare_refseq.py` — gene-call accuracy vs a reference GFF3
  (genes matched by contig+strand+stop; exact vs start-codon error; missed; extra).
- `scripts/download_data.py` (NCBI genome + RefSeq GFF3), `scripts/run_cam.py` (run + evaluate).
- `tests/test_cam.py` — pytest on a synthetic genome with planted genes (no internet needed).
- `docs/decisions.md`, `docs/problems.md`, `docs/meetings.md` — dissertation evidence logs;
  append to them when making design decisions or hitting problems.

Not built yet: PAM agent (functional annotation), CAM <-> PAM feedback loop, MCP tool server,
Swiss-Prot function evaluation.

## Tooling

The project is managed with [`uv`](https://docs.astral.sh/uv/) (v0.11.x), targeting
**Python 3.12** (pinned in `.python-version`; `pyproject.toml` requires `>=3.12`). The
`.venv/` is uv-managed — don't hand-edit it or create venvs with other tools.

```bash
uv run main.py --genome G.fna --ref R.gff3   # run CAM (same as scripts/run_cam.py)
uv run python -c "..."    # run arbitrary Python inside the managed env
uv add <package>          # add a runtime dependency (writes pyproject.toml + uv.lock)
uv add --dev <package>    # add a dev-only dependency
uv sync                   # install/refresh the env from pyproject.toml / uv.lock
uv lock                   # regenerate the lockfile
```

`uv.lock` is committed; keep it in sync with `pyproject.toml`.

### Tests

pytest is a dev dependency; config is in `pyproject.toml` (`pythonpath = ["src"]`).
`uv run pytest -q`, or a single test: `uv run pytest tests/test_cam.py::test_name`.

## Architecture

GeneWeaver v3 is a Master's dissertation project (supervisor: Kevin Bryson), a modern
successor to AGMIAL (Bryson et al., NAR 2006). Two agents, approved by the supervisor:
**CAM** (genome -> gene models -> proteins) and **PAM** (proteins -> name, EC, GO with
evidence), with PAM able to send suspect gene models back to CAM. Broader context: a
**multi-agent AI system for bioinformatics data integration** — an autonomous agent pipeline, *not* a web app. The
intended design is a set of cooperating agents that query bioinformatics databases
(UniProt, PDB, NCBI, Ensembl), coordinated by LLM orchestration across multi-step
analysis workflows. Expected stack: LangChain/LangGraph for orchestration, the A2A
(agent-to-agent) protocol and MCP for agent/tool communication, and ESM-2 protein
language models. Keep this direction in mind when structuring new code, but verify
against what actually exists before assuming any of it is implemented.
