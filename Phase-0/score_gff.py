"""Score one or more GFF3 annotations against the frozen E. coli K-12 MG1655 reference.

Usage (from Phase-0/):
    python3 score_gff.py --pred bakta=bakta/result.gff3 \
                         --pred pipeline=pipeline/pipeline_swissprot.gff3 \
                         --outdir results

Defaults point at the frozen reference in reference/. See README.md for the full rules.

Gene matching
    A predicted CDS matches a reference gene when it is on the same sequence and strand and has
    the same stop coordinate (end on +, start on -). A reference "gene" is one strand+stop locus;
    reference CDS that share a stop (alternative starts) form one gene.
    Reference pseudogenes are left out of recall. A prediction whose stop lies inside a reference
    pseudogene is not a false positive; it is counted separately.
    Predicted pseudogenes (pseudo=True or a pseudogene= attribute) are not counted as CDS predictions.

Metrics
    recall      matched reference genes / reference genes
    precision   matched predicted loci / predicted loci (excluding those on reference pseudogenes)
    start       matched genes whose predicted start equals a reference start
    function    for matched genes: gene symbol, EC (exact, full ECs only), normalised product,
                hypothetical-where-known, wrong function
    hypothetical share and pseudogene count over the whole prediction

No LLM is used anywhere: every comparison is a fixed string rule.
Uses only the Python standard library.
"""
import argparse
import csv
import random
import re
import urllib.parse
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_REFERENCE_GFF = HERE / "reference/ncbi_dataset/data/GCF_000005845.2/genomic.gff"
DEFAULT_REFERENCE_FUNCTION = HERE / "reference/reference_function.tsv"

SAMPLE_SIZE = 50
SAMPLE_SEED = 42

# A product is hypothetical when every word in it is one of these, a DUF/UPF/Pfam number,
# or a gene name such as YaaJ. "Uncharacterized transporter YaaJ" is not hypothetical: it names a function.
UNINFORMATIVE_WORDS = {
    "hypothetical", "uncharacterized", "uncharacterised", "putative", "probable", "predicted",
    "conserved", "protein", "domain", "containing", "family", "unknown", "function", "of",
}
UNINFORMATIVE_TOKEN = re.compile(r"(duf|upf|pf)\d+", re.IGNORECASE)
GENE_NAME_TOKEN = re.compile(r"[A-Za-z][a-z]{2}[A-Z]\d*")
# Words dropped before comparing product names
PRODUCT_STOPWORDS = {"putative", "probable", "possible", "predicted", "potential", "likely", "protein"}
FULL_EC = re.compile(r"^\d+\.\d+\.\d+\.n?\d+$")


# ---------------------------------------------------------------- GFF3 reading

def parse_attributes(text):
    """Parse a GFF3 attribute column; values are percent-decoded, multi-values kept as one string."""
    attributes = {}
    for item in text.strip().strip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key.strip()] = value
    return attributes


def decode(value):
    return urllib.parse.unquote(value or "")


def read_gff_features(path):
    """Yield (seqid, type, start, end, strand, attributes) for every feature line (stops at ##FASTA)."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            yield fields[0], fields[2], int(fields[3]), int(fields[4]), fields[6], parse_attributes(fields[8])


def is_pseudo(attributes):
    return attributes.get("pseudo", "").lower() == "true" or "pseudogene" in attributes


def ec_numbers(attributes):
    """All EC numbers on a feature: Dbxref=EC:..., ec_number=..., EC_number=..."""
    values = []
    for ref in attributes.get("Dbxref", "").split(","):
        ref = decode(ref)
        if ref.startswith("EC:"):
            values.append(ref[3:])
    for key in ("ec_number", "EC_number", "ec_numbers"):
        values += [decode(v) for v in re.split(r"[,;]", attributes.get(key, ""))]
    return sorted({value.strip() for value in values if value.strip()})


def stop_of(strand, start, end):
    return end if strand == "+" else start


def start_of(strand, start, end):
    return start if strand == "+" else end


# ---------------------------------------------------------------- reference

def load_reference(gff_path, function_path):
    """Return (genes, pseudo_spans).

    genes: {(seqid, strand, stop): gene} for non-pseudo reference CDS, one gene per stop locus.
    pseudo_spans: {(seqid, strand): [(start, end)]} for pseudogene features and pseudo CDS.
    """
    function = {}
    with open(function_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            function[row["cds_id"]] = row

    # Merge split CDS rows into one span per CDS ID
    cds = {}
    pseudo_spans = defaultdict(list)
    for seqid, ftype, start, end, strand, attributes in read_gff_features(gff_path):
        if ftype == "pseudogene":
            pseudo_spans[(seqid, strand)].append((start, end))
            continue
        if ftype != "CDS":
            continue
        record = cds.setdefault(attributes["ID"], {
            "seqid": seqid, "strand": strand, "start": start, "end": end, "pseudo": is_pseudo(attributes),
        })
        record["start"] = min(record["start"], start)
        record["end"] = max(record["end"], end)

    genes = {}
    for cds_id, record in cds.items():
        if record["pseudo"]:
            pseudo_spans[(record["seqid"], record["strand"])].append((record["start"], record["end"]))
            continue
        row = function.get(cds_id, {})
        key = (record["seqid"], record["strand"], stop_of(record["strand"], record["start"], record["end"]))
        gene = genes.setdefault(key, {
            "cds_ids": [], "starts": set(), "symbols": set(), "products": [], "ecs": set(), "locus_tags": set(),
        })
        gene["cds_ids"].append(cds_id)
        gene["starts"].add(start_of(record["strand"], record["start"], record["end"]))
        gene["locus_tags"].add(row.get("locus_tag", ""))
        for symbol in (row.get("refseq_gene"), row.get("uniprot_gene")):
            if symbol:
                gene["symbols"].add(symbol.lower())
        if row.get("refseq_product"):
            gene["products"].append(decode(row["refseq_product"]))
        if row.get("uniprot_protein_name"):
            gene["products"].append(row["uniprot_protein_name"].split(" (")[0])  # recommended name only
        gene["ecs"].update(ec for ec in row.get("ec_numbers", "").split(";") if ec)
    return genes, pseudo_spans


# ---------------------------------------------------------------- predictions

def load_predictions(path):
    """Return (loci, pseudo_count, duplicate_count).

    loci: {(seqid, strand, stop): prediction}, first CDS kept when several share a stop.
    CDS rows that share an ID (e.g. ribosomal slippage) are merged into one span first.
    """
    cds, pseudo_ids = {}, set()
    for number, (seqid, ftype, start, end, strand, attributes) in enumerate(read_gff_features(path)):
        feature_id = attributes.get("ID") or f"line{number}"
        if ftype == "pseudogene" or (ftype == "CDS" and is_pseudo(attributes)):
            pseudo_ids.add(feature_id)
            continue
        if ftype != "CDS":
            continue
        record = cds.setdefault(feature_id, {"seqid": seqid, "strand": strand, "start": start, "end": end,
                                             "attributes": attributes})
        record["start"] = min(record["start"], start)
        record["end"] = max(record["end"], end)

    loci, duplicates = {}, 0
    for feature_id, record in cds.items():
        strand, attributes = record["strand"], record["attributes"]
        key = (record["seqid"], strand, stop_of(strand, record["start"], record["end"]))
        if key in loci:
            duplicates += 1
            continue
        loci[key] = {
            "id": decode(attributes.get("ID", "")),
            "start": start_of(strand, record["start"], record["end"]),
            "product": decode(attributes.get("product", "")),
            "gene": decode(attributes.get("gene", "")),
            "ecs": ec_numbers(attributes),
        }
    return loci, len(pseudo_ids), duplicates


# ---------------------------------------------------------------- comparisons

def is_hypothetical(product):
    for word in re.split(r"[^A-Za-z0-9]+", product):
        if not word:
            continue
        if word.lower() in UNINFORMATIVE_WORDS or UNINFORMATIVE_TOKEN.fullmatch(word) or GENE_NAME_TOKEN.fullmatch(word):
            continue
        return False
    return True


def normalise_product(product):
    words = re.sub(r"[^a-z0-9.]+", " ", product.lower()).split()
    return " ".join(word for word in words if word not in PRODUCT_STOPWORDS)


def compare_function(prediction, gene):
    """Classify the function of one matched gene. Returns a dict of flags and a category."""
    full_ref_ecs = {ec for ec in gene["ecs"] if FULL_EC.match(ec)}
    full_pred_ecs = {ec for ec in prediction["ecs"] if FULL_EC.match(ec)}
    ref_products = {normalise_product(p) for p in gene["products"]} - {""}

    symbol_match = bool(prediction["gene"]) and prediction["gene"].lower() in gene["symbols"]
    ec_match = bool(full_ref_ecs & full_pred_ecs)
    product_match = normalise_product(prediction["product"]) in ref_products
    pred_hypothetical = is_hypothetical(prediction["product"])
    ref_known = bool(gene["products"]) and not is_hypothetical(gene["products"][0])  # RefSeq product decides

    if pred_hypothetical and ref_known:
        category = "hypothetical_where_known"
    elif pred_hypothetical:
        category = "both_hypothetical"
    elif not ref_known:
        category = "named_where_reference_hypothetical"
    elif symbol_match or ec_match or product_match:
        category = "agrees"
    else:
        category = "wrong_function"

    return {
        "symbol_match": symbol_match,
        "ec_scorable": bool(full_ref_ecs),
        "ec_partial_only": bool(gene["ecs"]) and not full_ref_ecs,
        "ec_match": ec_match,
        "ec_set_equal": bool(full_ref_ecs) and full_ref_ecs == full_pred_ecs,
        "product_match": product_match,
        "ref_known": ref_known,
        "category": category,
    }


def on_pseudogene(key, pseudo_spans):
    seqid, strand, stop = key
    return any(start <= stop <= end for start, end in pseudo_spans.get((seqid, strand), []))


# ---------------------------------------------------------------- scoring

def score(name, pred_path, genes, pseudo_spans, outdir):
    loci, pred_pseudo, duplicates = load_predictions(pred_path)

    matched_keys = [key for key in genes if key in loci]
    unmatched_pred = [key for key in loci if key not in genes]
    on_pseudo = [key for key in unmatched_pred if on_pseudogene(key, pseudo_spans)]
    false_pos = [key for key in unmatched_pred if not on_pseudogene(key, pseudo_spans)]

    gene_rows, comparisons = [], {}
    for key, gene in sorted(genes.items(), key=lambda item: item[0][2]):
        prediction = loci.get(key)
        row = {
            "seqid": key[0], "strand": key[1], "stop": key[2],
            "locus_tag": ";".join(sorted(t for t in gene["locus_tags"] if t)),
            "ref_start": ";".join(str(s) for s in sorted(gene["starts"])),
            "ref_symbols": ";".join(sorted(gene["symbols"])),
            "ref_product": gene["products"][0] if gene["products"] else "",
            "ref_ec": ";".join(sorted(gene["ecs"])),
            "status": "matched" if prediction else "missed",
            "pred_id": "", "pred_start": "", "start_exact": "",
            "pred_gene": "", "pred_product": "", "pred_ec": "",
            "symbol_match": "", "ec_match": "", "product_match": "", "category": "missed",
        }
        if prediction:
            result = compare_function(prediction, gene)
            comparisons[key] = result
            row.update({
                "pred_id": prediction["id"],
                "pred_start": prediction["start"],
                "start_exact": prediction["start"] in gene["starts"],
                "pred_gene": prediction["gene"],
                "pred_product": prediction["product"],
                "pred_ec": ";".join(prediction["ecs"]),
                "symbol_match": result["symbol_match"],
                "ec_match": result["ec_match"] if result["ec_scorable"] else "n/a",
                "product_match": result["product_match"],
                "category": result["category"],
            })
        gene_rows.append(row)

    write_tsv(outdir / f"{name}_genes.tsv", gene_rows)
    write_tsv(outdir / f"{name}_errors.tsv",
              [row for row in gene_rows if row["category"] in ("missed", "hypothetical_where_known")])
    write_tsv(outdir / f"{name}_false_positives.tsv", [
        {"seqid": k[0], "strand": k[1], "stop": k[2], "pred_id": loci[k]["id"],
         "pred_start": loci[k]["start"], "pred_product": loci[k]["product"]} for k in sorted(false_pos, key=lambda k: k[2])
    ])
    matched_rows = [row for row in gene_rows if row["status"] == "matched"]
    sample = random.Random(SAMPLE_SEED).sample(matched_rows, min(SAMPLE_SIZE, len(matched_rows)))
    write_tsv(outdir / f"{name}_handcheck_sample.tsv", [handcheck_row(n, row) for n, row in enumerate(sample, 1)])

    n_genes, n_matched = len(genes), len(matched_keys)
    n_pred_scored = len(loci) - len(on_pseudo)
    results = list(comparisons.values())
    count = lambda flag: sum(1 for r in results if r[flag])
    categories = defaultdict(int)
    for r in results:
        categories[r["category"]] += 1
    ec_scorable = count("ec_scorable")
    ref_known = count("ref_known")
    hypothetical_preds = sum(1 for p in loci.values() if is_hypothetical(p["product"]))

    return {
        "tool": name,
        "reference_genes": n_genes,
        "predicted_cds": len(loci),
        "matched": n_matched,
        "recall": n_matched / n_genes,
        "false_positives": len(false_pos),
        "on_reference_pseudogenes": len(on_pseudo),
        "precision": n_matched / n_pred_scored if n_pred_scored else 0.0,
        "start_exact": count_true(gene_rows, "start_exact") / n_matched if n_matched else 0.0,
        "symbol_match": count("symbol_match") / n_matched if n_matched else 0.0,
        "ec_scorable_genes": ec_scorable,
        "ec_match": count("ec_match") / ec_scorable if ec_scorable else 0.0,
        "ec_set_equal": count("ec_set_equal") / ec_scorable if ec_scorable else 0.0,
        "ec_partial_only_genes": count("ec_partial_only"),
        "product_match": count("product_match") / n_matched if n_matched else 0.0,
        "ref_known_genes": ref_known,
        "hypothetical_where_known": categories["hypothetical_where_known"] / ref_known if ref_known else 0.0,
        "wrong_function": categories["wrong_function"] / ref_known if ref_known else 0.0,
        "agrees": categories["agrees"] / ref_known if ref_known else 0.0,
        "hypothetical_share": hypothetical_preds / len(loci) if loci else 0.0,
        "pseudogenes": pred_pseudo,
        "duplicate_stops": duplicates,
    }


def handcheck_row(number, row):
    """One line of the hand-check sheet: reference and prediction side by side, then the scorer's calls."""
    return {
        "n": number, "locus_tag": row["locus_tag"], "strand": row["strand"], "stop": row["stop"],
        "ref_start": row["ref_start"], "pred_start": row["pred_start"],
        "ref_symbols": row["ref_symbols"], "pred_gene": row["pred_gene"],
        "ref_ec": row["ref_ec"], "pred_ec": row["pred_ec"],
        "ref_product": row["ref_product"], "pred_product": row["pred_product"],
        "start_exact": row["start_exact"], "symbol_match": row["symbol_match"], "ec_match": row["ec_match"],
        "product_match": row["product_match"], "category": row["category"],
        "hand_same_gene": "", "hand_same_function": "", "hand_category_ok": "", "hand_note": "",
    }


def count_true(rows, column):
    return sum(1 for row in rows if row[column] is True)


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as out:
        if not rows:
            out.write("")
            return
        writer = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------- report

TABLE_ROWS = [
    ("Reference genes", "reference_genes", "int"),
    ("Predicted CDS", "predicted_cds", "int"),
    ("Matched (same strand + stop)", "matched", "int"),
    ("Recall", "recall", "pct"),
    ("Precision", "precision", "pct"),
    ("False positives", "false_positives", "int"),
    ("Predictions on reference pseudogenes", "on_reference_pseudogenes", "int"),
    ("Start exact (of matched)", "start_exact", "pct"),
    ("Gene symbol match (of matched)", "symbol_match", "pct"),
    ("EC exact match (of matched with full ref EC)", "ec_match", "pct"),
    ("EC set identical (same denominator)", "ec_set_equal", "pct"),
    ("Matched genes with full ref EC", "ec_scorable_genes", "int"),
    ("Matched genes with only partial ref EC (not scored)", "ec_partial_only_genes", "int"),
    ("Product match, normalised (of matched)", "product_match", "pct"),
    ("Function agrees: symbol, EC or product (of matched, ref known)", "agrees", "pct"),
    ("Hypothetical where reference knows function", "hypothetical_where_known", "pct"),
    ("Wrong function (of matched, ref known)", "wrong_function", "pct"),
    ("Matched genes with known reference function", "ref_known_genes", "int"),
    ("Hypothetical share (of all predicted CDS)", "hypothetical_share", "pct"),
    ("Pseudogenes predicted", "pseudogenes", "int"),
]


def fmt(value, kind):
    return f"{100 * value:.1f}%" if kind == "pct" else f"{value:,}"


def write_report(metrics, outdir):
    write_tsv(outdir / "metrics.tsv", [{k: (round(v, 4) if isinstance(v, float) else v) for k, v in m.items()}
                                       for m in metrics])
    lines = ["| Metric | " + " | ".join(m["tool"] for m in metrics) + " |",
             "|---|" + "---:|" * len(metrics)]
    for label, key, kind in TABLE_ROWS:
        lines.append(f"| {label} | " + " | ".join(fmt(m[key], kind) for m in metrics) + " |")
    table = "\n".join(lines) + "\n"
    (outdir / "metrics.md").write_text(table)
    return table


def write_combined_errors(names, outdir):
    """One row per reference gene that any tool missed or called hypothetical, with every tool's outcome."""
    per_tool = {}
    for name in names:
        with open(outdir / f"{name}_genes.tsv") as fh:
            per_tool[name] = {(row["seqid"], row["strand"], row["stop"]): row for row in csv.DictReader(fh, delimiter="\t")}
    first = per_tool[names[0]]
    rows = []
    for key, ref in first.items():
        outcomes = {name: per_tool[name][key] for name in names}
        if not any(o["category"] in ("missed", "hypothetical_where_known") for o in outcomes.values()):
            continue
        row = {"seqid": key[0], "strand": key[1], "stop": key[2], "locus_tag": ref["locus_tag"],
               "ref_symbols": ref["ref_symbols"], "ref_product": ref["ref_product"], "ref_ec": ref["ref_ec"]}
        for name, outcome in outcomes.items():
            row[f"{name}_category"] = outcome["category"]
            row[f"{name}_product"] = outcome["pred_product"]
        row["tools_failing"] = sum(1 for o in outcomes.values() if o["category"] in ("missed", "hypothetical_where_known"))
        rows.append(row)
    rows.sort(key=lambda row: (-row["tools_failing"], int(row["stop"])))
    write_tsv(outdir / "errors_all_tools.tsv", rows)
    return len(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--pred", action="append", required=True, metavar="NAME=GFF3",
                        help="prediction to score; repeat for several tools")
    parser.add_argument("--reference-gff", type=Path, default=DEFAULT_REFERENCE_GFF)
    parser.add_argument("--reference-function", type=Path, default=DEFAULT_REFERENCE_FUNCTION)
    parser.add_argument("--outdir", type=Path, default=HERE / "results")
    args = parser.parse_args()

    genes, pseudo_spans = load_reference(args.reference_gff, args.reference_function)
    metrics = []
    for item in args.pred:
        name, path = item.split("=", 1)
        metrics.append(score(name, Path(path), genes, pseudo_spans, args.outdir))
    print(write_report(metrics, args.outdir))
    if len(metrics) > 1:
        count = write_combined_errors([m["tool"] for m in metrics], args.outdir)
        print(f"{count} reference genes missed or called hypothetical by at least one tool -> errors_all_tools.tsv")


if __name__ == "__main__":
    main()
