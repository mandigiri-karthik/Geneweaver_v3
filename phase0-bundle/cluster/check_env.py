#!/usr/bin/env python3
"""Verify the Phase 0 environment before running anything.

Checks the Python packages, the DIAMOND and Bakta binaries, the repository layout
and every reference database (presence, size, and optionally md5 against
../manifest.tsv). Prints one table and exits non-zero if anything required is
missing, so a notebook cell fails loudly rather than halfway through a run.

    python3 cluster/check_env.py                 # fast: presence + size
    python3 cluster/check_env.py --verify-md5    # slow: full md5 (~2.8 GB)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REQUIRED_DBS = [
    "uniprot_sprot.fasta.gz",
    "uniprot_sprot_no_ecoli.fasta.gz",
    "Pfam-A.hmm",
    "Pfam-A.clans.tsv.gz",
    "pfam2go",
    "swissprot.dmnd",
    "swissprot_no_ecoli.dmnd",
]
REPO_PATHS = [
    "Structural_annotation/predict_genes.py",
    "Structural_annotation/homology_annotator.py",
    "Structural_annotation/pfam_annotator.py",
    "Structural_annotation/go_mapper.py",
    "Structural_annotation/annotation_combiner.py",
    "Phase-0/score_gff.py",
    "Phase-0/sanity_checks.py",
    "Phase-0/export_pipeline_gff3.py",
    "Phase-0/bakta/replicons.tsv",
    "Phase-0/reference/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna",
]

rows: list[tuple[str, str, str]] = []
problems: list[str] = []


def add(name: str, ok: bool | None, detail: str, required: bool = True) -> None:
    mark = "ok " if ok else ("-- " if ok is None else "FAIL")
    rows.append((mark, name, detail))
    if ok is False and required:
        problems.append(f"{name}: {detail}")


def md5(path: Path, chunk: int = 1 << 22) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def manifest() -> dict[str, str]:
    out: dict[str, str] = {}
    f = HERE.parent / "manifest.tsv"
    if not f.is_file():
        return out
    for line in f.read_text().splitlines():
        if line.startswith("#") or not line.strip() or line.startswith("file\t"):
            continue
        parts = line.split("\t")
        if len(parts) >= 3 and parts[2].strip() == "yes":
            out[parts[0]] = parts[1]
    return out


def human(n: int) -> str:
    for unit in ("B", "K", "M", "G"):
        if n < 1024 or unit == "G":
            return f"{n:.0f}{unit}" if unit == "B" else f"{n/1:.1f}{unit}" if False else f"{n:.1f}{unit}"
        n /= 1024.0
    return f"{n:.1f}G"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--paths", type=Path, default=HERE / "paths.json")
    p.add_argument("--verify-md5", action="store_true",
                   help="md5 every database (slow; ~2.8 GB)")
    a = p.parse_args()

    print(f"python {sys.version.split()[0]}  ({sys.executable})")

    if not a.paths.is_file():
        print(f"\nERROR: {a.paths} not found.\n"
              f"Copy cluster/paths.example.json to cluster/paths.json and edit it.",
              file=sys.stderr)
        return 2
    cfg = json.loads(a.paths.read_text())
    repo = Path(cfg["repo"]).expanduser()
    dbs = Path(cfg["databases"]).expanduser()
    bakta_db = Path(cfg.get("bakta_db", "")).expanduser() if cfg.get("bakta_db") else None

    # --- python packages -----------------------------------------------------
    for mod, attr in [("pyrodigal", "__version__"), ("pyhmmer", "__version__"), ("Bio", "__version__")]:
        try:
            m = __import__(mod)
            add(f"python: {mod}", True, getattr(m, attr, "?"))
        except ImportError as e:
            add(f"python: {mod}", False, f"not importable ({e})")

    # --- binaries ------------------------------------------------------------
    dpath = shutil.which("diamond")
    if dpath:
        v = subprocess.run([dpath, "--version"], capture_output=True, text=True).stdout.strip()
        add("binary: diamond", True, f"{v}  [{dpath}]")
    else:
        add("binary: diamond", False, "not on PATH - run cluster/bootstrap.sh and export PATH")

    bakta_env = cfg.get("bakta_env")
    if bakta_env and Path(bakta_env).is_dir():
        r = subprocess.run(["conda", "run", "-p", bakta_env, "bakta", "--version"],
                           capture_output=True, text=True)
        ok = r.returncode == 0
        add("binary: bakta", ok or None,
            (r.stdout.strip() or r.stderr.strip()[:80]) + f"  [-p {bakta_env}]", required=False)
    else:
        bpath = shutil.which("bakta")
        if bpath:
            v = subprocess.run([bpath, "--version"], capture_output=True, text=True).stdout.strip()
            add("binary: bakta", True, f"{v}  [{bpath}]", required=False)
        else:
            add("binary: bakta", None,
                "set \"bakta_env\" in paths.json (the env has a path, not a name)", required=False)

    # --- repository ----------------------------------------------------------
    add("repo root", repo.is_dir(), str(repo))
    for rel in REPO_PATHS:
        f = repo / rel
        note = ""
        if not f.is_file() and "Structural_annotation" in rel and rel.endswith(
                ("homology_annotator.py", "pfam_annotator.py", "go_mapper.py", "annotation_combiner.py")):
            note = "  <- on branch 'dev'; merge it into main"
        add(f"repo: {rel}", f.is_file(), (("" if f.is_file() else "missing") + note).strip() or "present")

    # --- databases -----------------------------------------------------------
    add("database dir", dbs.is_dir(), str(dbs))
    want = manifest() if a.verify_md5 else {}
    for name in REQUIRED_DBS:
        f = dbs / name
        if not f.is_file():
            add(f"db: {name}", False, "missing - run scripts/fetch_databases.sh")
            continue
        detail = f"{f.stat().st_size / 1024**2:.0f} MB"
        if a.verify_md5 and name in want:
            got = md5(f)
            if got == want[name]:
                detail += "  md5 ok"
            else:
                add(f"db: {name}", False, f"md5 mismatch (expected {want[name]}, got {got})")
                continue
        add(f"db: {name}", True, detail)

    if bakta_db is not None:
        v = bakta_db / "version.json"
        if v.is_file():
            try:
                info = json.loads(v.read_text())
                add("bakta db", True, f"v{info.get('major')}.{info.get('minor')} {info.get('type','')}".strip(),
                    required=False)
            except Exception:
                add("bakta db", True, str(bakta_db), required=False)
        else:
            add("bakta db", bakta_db.is_dir(), str(bakta_db), required=False)

    # --- storage -------------------------------------------------------------
    # The IDA pod mounts five PVCs; all but /mnt/primary are 95-100% full, so writing
    # anywhere else fails late and confusingly.
    for target in {repo, dbs}:
        try:
            st = shutil.disk_usage(target if target.exists() else target.parent)
            free_gb = st.free / 1024**3
            pct = 100 * st.used / st.total
            add(f"disk: {target}", free_gb > 20, f"{free_gb:.0f} GB free ({pct:.0f}% used)",
                required=False)
        except OSError as e:
            add(f"disk: {target}", None, str(e), required=False)

    # --- report --------------------------------------------------------------
    w = max(len(r[1]) for r in rows)
    print()
    for mark, name, detail in rows:
        print(f"  [{mark}] {name:<{w}}  {detail}")

    if problems:
        print(f"\n{len(problems)} problem(s) must be fixed before running Phase 0:", file=sys.stderr)
        for s in problems:
            print(f"  - {s}", file=sys.stderr)
        return 1
    print("\nEnvironment OK.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
