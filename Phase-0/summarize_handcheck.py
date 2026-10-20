"""Summarise filled-in hand-check sheets: how often the scorer's rules agree with a human judgement.

Usage (from Phase-0/):
    python3 summarize_handcheck.py results/*_handcheck_sample.tsv

Reads the three hand columns (Y or N; blank = not checked yet):
    hand_same_gene      the predicted CDS really is the reference gene
    hand_same_function  the predicted function really is the reference function
    hand_category_ok    the scorer's category is right

Uses only the Python standard library.
"""
import csv
import sys
from pathlib import Path

HAND_COLUMNS = ["hand_same_gene", "hand_same_function", "hand_category_ok"]


def summarise(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    print(f"\n{Path(path).name}: {len(rows)} cases")

    for column in HAND_COLUMNS:
        answers = [row[column].strip().upper()[:1] for row in rows if row[column].strip()]
        yes = answers.count("Y")
        share = f"{100 * yes / len(answers):.0f}%" if answers else "n/a"
        print(f"  {column:<20} checked {len(answers):>3}  yes {yes:>3}  ({share})")

    # Where the strict rules said "no function match" but the human says the function is the same
    too_strict = [row for row in rows if row["category"] == "wrong_function"
                  and row["hand_same_function"].strip().upper().startswith("Y")]
    flagged_wrong = sum(1 for row in rows if row["category"] == "wrong_function")
    print(f"  scorer said wrong_function: {flagged_wrong}; human says function is the same: {len(too_strict)}")
    for row in too_strict:
        print(f"    {row['locus_tag']:<8} ref '{row['ref_product']}'  vs  pred '{row['pred_product']}'")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for path in sys.argv[1:]:
        summarise(path)


if __name__ == "__main__":
    main()
