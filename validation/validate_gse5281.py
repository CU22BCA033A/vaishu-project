"""
validate_gse5281.py
-------------------
Tests the 7-gene panel on a dataset it has never seen: GSE5281.

Why: the panel was chosen using GSE48350. Any score measured on GSE48350
itself is too optimistic. GSE5281 comes from a different study and lab,
so it gives an honest check.

What this script does:
  1. Downloads the GSE5281 "series matrix" file from NCBI GEO (once).
  2. Keeps only the hippocampus samples (control and Alzheimer's).
  3. Picks one probe per panel gene: the probe with the highest average
     expression (the same rule used to build the panel).
  4. Scores every sample WITHOUT looking at the labels:
        - each gene is standardized within this dataset (z-score)
        - panel score = average of sign * z
  5. Only then compares the scores with the real labels:
        - AUC (0.5 = coin flip, 1.0 = perfect separation)
        - accuracy using a fixed rule decided in advance: score > 0
        - a permutation test (how often random labels do as well)
        - a bootstrap 95% interval for the AUC
  6. Saves the results in validation/results/.

Nothing about the panel (genes, directions) is changed after seeing the data.

Run from the project folder:
    python validation/validate_gse5281.py
"""

import csv
import gzip
import json
import math
import os
import random
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_FOLDER = os.path.join(HERE, "data")
RESULTS_FOLDER = os.path.join(HERE, "results")
MATRIX_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5281/matrix/"
              "GSE5281_series_matrix.txt.gz")
MATRIX_FILE = os.path.join(DATA_FOLDER, "GSE5281_series_matrix.txt.gz")
PANEL_FILE = os.path.join(HERE, "..", "static", "data", "panel_reference.csv")

# Probes on the Affymetrix HG-U133 Plus 2.0 chip (GPL570) for each panel gene.
# Found in the GPL570 annotation table on GEO. The two LINC genes have newer
# names than that table, so they were matched by their NCBI Gene ID:
#   LINC03007 = Gene ID 646588    = "LOC646588"    in the chip table
#   LINC02987 = Gene ID 101927151 = "LOC101927151" in the chip table
PROBES = {
    "ANKIB1":    ["224682_at", "224687_at", "243181_at"],
    "RAE1":      ["201558_at", "211318_s_at", "211319_at"],
    "SLC25A46":  ["212833_at", "226831_at"],
    "LINC03007": ["236977_at"],
    "ZNF621":    ["1558620_at", "230164_at", "241499_at"],
    "ONECUT3":   ["234109_x_at"],
    "LINC02987": ["1569478_s_at"],
}

random.seed(42)  # so the permutation test and bootstrap give the same numbers every run


# ---------------------------------------------------------------------------
# Small maths helpers
# ---------------------------------------------------------------------------

def mean(values):
    return sum(values) / len(values)


def sd(values):
    m = mean(values)
    return math.sqrt(sum((v - m) ** 2 for v in values) / (len(values) - 1))


def auc(scores, labels):
    """
    Area under the ROC curve = the chance that a random Alzheimer's sample
    has a higher score than a random control sample (ties count half).
    labels: 1 = Alzheimer's, 0 = control
    """
    ad = [s for s, y in zip(scores, labels) if y == 1]
    ctl = [s for s, y in zip(scores, labels) if y == 0]
    wins = 0.0
    for a in ad:
        for c in ctl:
            if a > c:
                wins += 1
            elif a == c:
                wins += 0.5
    return wins / (len(ad) * len(ctl))


# ---------------------------------------------------------------------------
# Step 1 and 2: download and read the data
# ---------------------------------------------------------------------------

def download_if_needed():
    if os.path.exists(MATRIX_FILE):
        return
    os.makedirs(DATA_FOLDER, exist_ok=True)
    print("Downloading GSE5281 (about 36 MB) ...")
    urllib.request.urlretrieve(MATRIX_URL, MATRIX_FILE)


def read_series_matrix():
    """
    Returns:
        samples: list of dicts {id, title, region, disease, cell_type, sex, age}
        values:  dict probe -> list of numbers (same order as samples)
    """
    wanted = set()
    for probe_list in PROBES.values():
        wanted.update(probe_list)

    info = {}            # header lines we care about
    characteristics = [] # all "!Sample_characteristics_ch1" lines
    values = {}
    in_table = False

    with gzip.open(MATRIX_FILE, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            parts = [p.strip().strip('"') for p in line.split("\t")]
            if in_table:
                probe = parts[0]
                if probe in wanted:
                    values[probe] = [float(x) if x not in ("", "null") else None for x in parts[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                characteristics.append(parts[1:])
            elif line.startswith("!Sample_title") or line.startswith("!Sample_geo_accession"):
                info[parts[0]] = parts[1:]

    # Turn "Organ Region: hippocampus" style lines into per-sample fields
    n = len(info["!Sample_title"])
    fields = [{} for _ in range(n)]
    for row in characteristics:
        for i, cell in enumerate(row):
            if ":" in cell:
                key, value = cell.split(":", 1)
                fields[i][key.strip().lower()] = value.strip()

    samples = []
    for i in range(n):
        samples.append({
            "id": info["!Sample_geo_accession"][i],
            "title": info["!Sample_title"][i],
            "region": fields[i].get("organ region", ""),
            "disease": fields[i].get("disease state", ""),
            "cell_type": fields[i].get("cell type", ""),
            "sex": fields[i].get("sex", ""),
            "age": fields[i].get("age", ""),
        })
    return samples, values


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def main():
    download_if_needed()
    samples, values = read_series_matrix()

    # Keep hippocampus samples only
    keep = [i for i, s in enumerate(samples) if "hippocampus" in s["region"].lower()]
    hip = [samples[i] for i in keep]
    labels = [0 if "normal" in s["disease"].lower() else 1 for s in hip]
    print("Hippocampus samples:", len(hip), "| control:", labels.count(0), "| Alzheimer's:", labels.count(1))
    print("Cell type:", sorted(set(s["cell_type"] for s in hip)))

    # MAS5 values are on a linear scale -> log2 (add 1 so 0 stays finite)
    def log2_row(probe):
        return [math.log2(max(values[probe][i], 0) + 1) for i in keep]

    # Read the panel (genes + direction). The control means/SDs from GSE48350
    # are NOT used here, because GSE5281 was processed differently (MAS5,
    # microdissected neurons): its numbers are on another scale.
    panel = []
    with open(PANEL_FILE) as f:
        for row in csv.DictReader(f):
            panel.append({"gene": row["gene"], "sign": int(row["sign"])})

    # Step 3: one probe per gene = highest average expression (label-free)
    chosen = {}
    gene_values = {}
    for p in panel:
        best_probe = None
        best_mean = None
        for probe in PROBES[p["gene"]]:
            if probe not in values:
                continue
            m = mean(log2_row(probe))
            if best_mean is None or m > best_mean:
                best_probe, best_mean = probe, m
        chosen[p["gene"]] = best_probe
        gene_values[p["gene"]] = log2_row(best_probe)

    # Step 4: within-dataset z-scores and panel score (labels NOT used)
    z = {}
    for p in panel:
        vals = gene_values[p["gene"]]
        m, s = mean(vals), sd(vals)
        z[p["gene"]] = [(v - m) / s for v in vals]

    scores = []
    for i in range(len(hip)):
        contributions = [p["sign"] * z[p["gene"]][i] for p in panel]
        scores.append(mean(contributions))

    # Step 5: compare with the true labels
    panel_auc = auc(scores, labels)
    predicted = [1 if s > 0 else 0 for s in scores]
    correct = sum(1 for a, b in zip(predicted, labels) if a == b)
    tp = sum(1 for a, b in zip(predicted, labels) if a == 1 and b == 1)
    tn = sum(1 for a, b in zip(predicted, labels) if a == 0 and b == 0)

    # Permutation test: shuffle labels 10,000 times
    better = 0
    shuffled = labels[:]
    for _ in range(10000):
        random.shuffle(shuffled)
        if auc(scores, shuffled) >= panel_auc:
            better += 1
    p_value = (better + 1) / (10000 + 1)

    # Bootstrap 95% interval for the AUC
    boot = []
    n = len(scores)
    while len(boot) < 2000:
        idx = [random.randrange(n) for _ in range(n)]
        ys = [labels[i] for i in idx]
        if 0 < sum(ys) < n:  # need both groups in the resample
            boot.append(auc([scores[i] for i in idx], ys))
    boot.sort()
    ci = (boot[int(0.025 * len(boot))], boot[int(0.975 * len(boot)) - 1])

    # Per-gene check: does each gene move in the expected direction?
    genes = []
    for p in panel:
        vals = gene_values[p["gene"]]
        ad_mean = mean([v for v, y in zip(vals, labels) if y == 1])
        ctl_mean = mean([v for v, y in zip(vals, labels) if y == 0])
        diff = ad_mean - ctl_mean
        gene_auc = auc([p["sign"] * v for v in vals], labels)
        genes.append({
            "gene": p["gene"],
            "probe": chosen[p["gene"]],
            "expected": "higher" if p["sign"] > 0 else "lower",
            "observed": "higher" if diff > 0 else "lower",
            "agrees": (diff > 0) == (p["sign"] > 0),
            "ad_minus_control_log2": round(diff, 3),
            "auc_in_expected_direction": round(gene_auc, 3),
            "mean_log2": round(mean(vals), 2),
        })

    # Sensitivity check: does a different probe choice change the answer?
    # (Reported only. The main result uses the rule fixed in advance.)
    def auc_with_probe_rule(rule):
        rule_scores = [0.0] * len(hip)
        for p in panel:
            rows = [log2_row(probe) for probe in PROBES[p["gene"]] if probe in values]
            if rule == "average of all probes":
                gene_row = [mean(column) for column in zip(*rows)]
            else:  # "most variable probe"
                gene_row = max(rows, key=sd)
            m, s = mean(gene_row), sd(gene_row)
            for i, v in enumerate(gene_row):
                rule_scores[i] += p["sign"] * (v - m) / s / len(panel)
        return round(auc(rule_scores, labels), 3)

    sensitivity_check = {
        "highest-expressed probe (main rule)": round(panel_auc, 3),
        "average of all probes": auc_with_probe_rule("average of all probes"),
        "most variable probe": auc_with_probe_rule("most variable probe"),
    }

    summary = {
        "dataset": "GSE5281",
        "tissue": "hippocampus (CA1 pyramidal neurons, laser-capture microdissected)",
        "platform": "Affymetrix HG-U133 Plus 2.0 (GPL570), MAS5 signal, log2(x+1)",
        "n_control": labels.count(0),
        "n_alzheimers": labels.count(1),
        "scoring": "within-dataset z-score per gene, panel score = mean(sign * z); labels not used",
        "auc": round(panel_auc, 3),
        "auc_ci95": [round(ci[0], 3), round(ci[1], 3)],
        "permutation_p": round(p_value, 4),
        "rule": "score > 0 = Alzheimer's-like (fixed in advance)",
        "accuracy": round(correct / n, 3),
        "correct": correct,
        "total": n,
        "sensitivity": round(tp / labels.count(1), 3),
        "specificity": round(tn / labels.count(0), 3),
        "genes_in_expected_direction": sum(1 for g in genes if g["agrees"]),
        "genes": genes,
        "probe_choice_check_auc": sensitivity_check,
    }

    os.makedirs(RESULTS_FOLDER, exist_ok=True)
    with open(os.path.join(RESULTS_FOLDER, "gse5281_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    with open(os.path.join(RESULTS_FOLDER, "gse5281_scores.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample", "title", "group", "panel_score"])
        for s, y, sc in zip(hip, labels, scores):
            w.writerow([s["id"], s["title"], "Alzheimer" if y else "control", round(sc, 4)])

    # The 7 panel genes for every sample, in the layout the website accepts.
    # Students can upload this file in "batch mode" to see the same result.
    with open(os.path.join(RESULTS_FOLDER, "gse5281_panel_genes.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample", "group"] + [p["gene"] for p in panel])
        for i, (s, y) in enumerate(zip(hip, labels)):
            w.writerow([s["id"], "Alzheimer" if y else "control"]
                       + [round(gene_values[p["gene"]][i], 4) for p in panel])

    # Print a short report
    print()
    print("AUC: %.3f  (95%% CI %.3f - %.3f), permutation p = %.4f" % (panel_auc, ci[0], ci[1], p_value))
    print("Accuracy with score > 0: %d/%d = %.0f%%  (sensitivity %.0f%%, specificity %.0f%%)" % (
        correct, n, 100 * correct / n, 100 * summary["sensitivity"], 100 * summary["specificity"]))
    print("Genes moving in the expected direction: %d of %d" % (summary["genes_in_expected_direction"], len(genes)))
    print("AUC with other probe choices:", sensitivity_check)
    for g in genes:
        print("  %-10s %-13s expected %-6s observed %-6s diff %+.3f  gene AUC %.2f  mean log2 %.1f" % (
            g["gene"], g["probe"], g["expected"], g["observed"], g["ad_minus_control_log2"],
            g["auc_in_expected_direction"], g["mean_log2"]))


if __name__ == "__main__":
    main()
