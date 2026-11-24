"""
analysis.py
-----------
All the "science" of the project lives in this file.
It does NOT know anything about websites - it only works with text and numbers.
That makes it easy to read and easy to test on its own.

The idea in 4 steps:
  1. Read the panel reference (7 genes, their direction, and the mean / SD
     measured in healthy control brains).
  2. Read the user's gene expression table and pick out only the panel genes.
  3. For every gene, compute a z-score:  z = (value - mean) / sd
     The mean and SD come either from the GSE48350 controls ("reference"
     mode) or from the uploaded samples themselves ("batch" mode).
     Multiply z by the gene's sign (+1 = higher in Alzheimer's,
     -1 = lower in Alzheimer's). This is the gene's "contribution".
  4. The panel score is the average contribution. If it is above the cutoff,
     the sample looks more like the Alzheimer's group than the control group.

IMPORTANT: this is an educational project, NOT a medical test.
"""

import csv
import io
import math
import os

# Folder where the default data files are stored.
DATA_FOLDER = os.path.join(os.path.dirname(__file__), "static", "data")
PANEL_FILE = os.path.join(DATA_FOLDER, "panel_reference.csv")
EXAMPLE_FILE = os.path.join(DATA_FOLDER, "example_expression.csv")

# 95th percentile of the 18 reference control scores from GSE48350.
DEFAULT_CUTOFF = 0.46

# A sample needs at least 80% of the panel genes to get a score.
MIN_GENE_FRACTION = 0.8


class AnalysisError(Exception):
    """A friendly error message that we can show to the user."""
    pass


# ---------------------------------------------------------------------------
# Small helper functions
# ---------------------------------------------------------------------------

def is_number(text):
    """Return True if the text can be turned into a number, e.g. '2.5'."""
    try:
        value = float(text)
    except (ValueError, TypeError):
        return False
    # float("nan") and float("inf") work, but they are not useful numbers
    return not math.isnan(value) and not math.isinf(value)


def read_file_text(path):
    """Read a whole text file and return it as one string."""
    with open(path, encoding="utf-8") as f:
        return f.read()


def read_table(text):
    """
    Turn CSV / TSV text into a list of rows (each row is a list of strings).
    We guess the separator: tab, semicolon or comma.
    """
    lines = [line for line in text.splitlines() if line.strip() != ""]
    if len(lines) == 0:
        return []

    first_line = lines[0]
    if "\t" in first_line:
        separator = "\t"
    elif first_line.count(";") > first_line.count(","):
        separator = ";"
    else:
        separator = ","

    reader = csv.reader(io.StringIO("\n".join(lines)), delimiter=separator)
    rows = []
    for row in reader:
        rows.append([cell.strip() for cell in row])
    return rows


def percentile(values, p):
    """
    The p-th percentile of a list of numbers (p between 0 and 1).
    Example: percentile([1, 2, 3, 4, 5], 0.5) -> 3
    """
    numbers = sorted(values)
    position = (len(numbers) - 1) * p
    low = math.floor(position)
    high = math.ceil(position)
    fraction = position - low
    return numbers[low] + (numbers[high] - numbers[low]) * fraction


def average(values):
    return sum(values) / len(values)


# ---------------------------------------------------------------------------
# Step 1: the panel reference
# ---------------------------------------------------------------------------

def parse_panel(text):
    """
    Read the panel reference. Expected columns:
        gene, sign, control_mean, control_sd
    Returns a list of dictionaries, one per gene.
    """
    rows = read_table(text)

    # Skip the header row if the 3rd column is not a number
    if len(rows) > 0 and len(rows[0]) > 2 and not is_number(rows[0][2]):
        rows = rows[1:]

    panel = []
    for row in rows:
        if len(row) < 4:
            continue
        gene = row[0].upper()
        sign_text = row[1].lower()
        if sign_text.startswith("-") or "down" in sign_text or "lower" in sign_text:
            sign = -1
        else:
            sign = 1

        if not is_number(row[2]) or not is_number(row[3]):
            continue
        mean = float(row[2])
        sd = float(row[3])
        if gene == "" or sd <= 0:
            continue

        panel.append({"gene": gene, "sign": sign, "mean": mean, "sd": sd})

    if len(panel) < 2:
        raise AnalysisError(
            "The panel reference needs at least 2 genes in the format "
            "gene,sign,control_mean,control_sd."
        )
    return panel


def default_panel_text():
    """The 7-gene panel that ships with the project."""
    return read_file_text(PANEL_FILE)


# ---------------------------------------------------------------------------
# Optional: reference scores (used to calculate the cutoff)
# ---------------------------------------------------------------------------

def parse_reference_scores(text):
    """
    Read panel_scores.csv (columns: sample, group, panel_score).
    Returns a dictionary with the control scores, the Alzheimer's scores
    and the cutoff (95th percentile of the control scores).
    """
    rows = read_table(text)
    if len(rows) < 2:
        raise AnalysisError("The reference scores file needs a header row and some data rows.")

    header = [name.lower() for name in rows[0]]

    # Find the right columns (fall back to sensible defaults)
    score_col = header.index("panel_score") if "panel_score" in header else len(header) - 1
    group_col = header.index("group") if "group" in header else 1

    controls = []
    alzheimers = []
    for row in rows[1:]:
        if len(row) <= max(score_col, group_col):
            continue
        if not is_number(row[score_col]):
            continue
        score = float(row[score_col])
        if "control" in row[group_col].lower():
            controls.append(score)
        else:
            alzheimers.append(score)

    if len(controls) < 5:
        raise AnalysisError(
            "The reference scores need at least 5 control samples "
            "(group column containing the word 'control') to calculate a cutoff."
        )

    return {
        "controls": controls,
        "alzheimers": alzheimers,
        "cutoff": percentile(controls, 0.95),
    }


# ---------------------------------------------------------------------------
# Step 2: the user's expression data
# ---------------------------------------------------------------------------

GENE_COLUMN_NAMES = ["gene", "genes", "symbol", "gene_symbol", "gene.symbol", "genesymbol"]


def parse_expression(text, panel):
    """
    Read the expression table. Two layouts are allowed:

      A) genes as rows      gene,sample1,sample2,...
                            APP,5.1,5.3,...

      B) samples as rows    sample,group,ANKIB1,RAE1,...
                            S1,control,0.4,2.6,...

    Returns a list of samples. Each sample is a dictionary:
        {"name": "sample1", "group": "", "values": {"ANKIB1": 0.4, ...}}
    """
    rows = read_table(text)
    if len(rows) < 2:
        raise AnalysisError("The expression data needs a header row and at least one more row.")

    panel_genes = [p["gene"] for p in panel]
    header = rows[0]

    # Which column holds the gene names? (layout A)
    gene_col = 0
    for i, name in enumerate(header):
        if name.lower() in GENE_COLUMN_NAMES:
            gene_col = i
            break

    # Count how many panel genes appear in the header vs. in the gene column
    genes_in_header = 0
    for name in header:
        if name.upper() in panel_genes:
            genes_in_header += 1

    genes_in_column = 0
    for row in rows[1:]:
        if len(row) > gene_col and row[gene_col].upper() in panel_genes:
            genes_in_column += 1

    if genes_in_header == 0 and genes_in_column == 0:
        raise AnalysisError(
            "None of the panel genes were found. Check that the gene symbols "
            "(e.g. " + ", ".join(panel_genes[:3]) + ") are in your file."
        )

    if genes_in_header > genes_in_column:
        return parse_samples_as_rows(rows, panel_genes)
    return parse_genes_as_rows(rows, panel_genes, gene_col)


def parse_samples_as_rows(rows, panel_genes):
    """Layout B: one row per sample, one column per gene."""
    header = rows[0]

    name_col = 0
    group_col = None
    for i, name in enumerate(header):
        if name.lower() in ("sample", "id", "sample_id"):
            name_col = i
        if name.lower() == "group":
            group_col = i

    samples = []
    for row in rows[1:]:
        values = {}
        for i, column_name in enumerate(header):
            gene = column_name.upper()
            if gene in panel_genes and i < len(row) and is_number(row[i]):
                values[gene] = float(row[i])

        group = ""
        if group_col is not None and group_col < len(row):
            group = row[group_col]

        samples.append({"name": row[name_col], "group": group, "values": values})
    return samples


def parse_genes_as_rows(rows, panel_genes, gene_col):
    """
    Layout A: one row per gene (or per probe), one column per sample.
    If a gene appears on several rows (several probes), we keep the row
    with the highest average expression - the same rule that was used
    to build the reference.
    """
    header = rows[0]
    data_rows = rows[1:]

    # Sample columns = every column (except the gene column) that holds numbers
    sample_cols = []
    for i in range(len(header)):
        if i == gene_col:
            continue
        numeric = 0
        checked = 0
        for row in data_rows[:20]:
            if i < len(row):
                checked += 1
                if is_number(row[i]):
                    numeric += 1
        if checked > 0 and numeric >= checked / 2:
            sample_cols.append(i)

    if len(sample_cols) == 0:
        raise AnalysisError("No columns with numbers were found in the expression data.")

    # If the header itself is numbers, the file has no header -> make names
    if is_number(header[sample_cols[0]]):
        names = ["sample" + str(n + 1) for n in range(len(sample_cols))]
        data_rows = rows
    else:
        names = [header[i] for i in sample_cols]

    # Optional row called "group" (e.g. group,control,AD,AD,...)
    groups = ["" for _ in sample_cols]
    for row in data_rows:
        if len(row) > gene_col and row[gene_col].lower() == "group":
            groups = [row[i] if i < len(row) else "" for i in sample_cols]

    # best_row[gene] = (average, row) for the strongest probe of that gene
    best_row = {}
    for row in data_rows:
        if len(row) <= gene_col:
            continue
        gene = row[gene_col].upper()
        if gene not in panel_genes:
            continue
        numbers = [float(row[i]) for i in sample_cols if i < len(row) and is_number(row[i])]
        if len(numbers) == 0:
            continue
        row_average = average(numbers)
        if gene not in best_row or row_average > best_row[gene][0]:
            best_row[gene] = (row_average, row)

    samples = []
    for n, col in enumerate(sample_cols):
        values = {}
        for gene, (row_average, row) in best_row.items():
            if col < len(row) and is_number(row[col]):
                values[gene] = float(row[col])
        samples.append({"name": names[n], "group": groups[n], "values": values})
    return samples


# ---------------------------------------------------------------------------
# Step 3 and 4: scoring
# ---------------------------------------------------------------------------
#
# Two ways to decide what "normal" looks like for each gene:
#
#   "reference" mode: use the control mean and SD measured in GSE48350.
#                     Only fair if your data was processed exactly like it.
#
#   "batch" mode:     use the mean and SD of YOUR uploaded samples.
#                     This removes most lab/processing differences, but it
#                     needs a batch of at least 10 samples that contains
#                     both kinds of samples. Labels are NOT used for this.

MIN_BATCH_SIZE = 10


def reference_centers(panel):
    """Mean and SD for each gene, taken from the panel reference file."""
    centers = {}
    for p in panel:
        centers[p["gene"]] = (p["mean"], p["sd"])
    return centers


def standard_deviation(values):
    m = average(values)
    return math.sqrt(sum((v - m) ** 2 for v in values) / (len(values) - 1))


def batch_centers(panel, samples):
    """Mean and SD for each gene, calculated from the uploaded samples."""
    centers = {}
    for p in panel:
        values = [s["values"][p["gene"]] for s in samples if p["gene"] in s["values"]]
        if len(values) < 3:
            continue
        sd = standard_deviation(values)
        if sd > 0:
            centers[p["gene"]] = (average(values), sd)
    return centers


def score_sample(sample, panel, centers):
    """
    Compare one sample with the panel.
    centers[gene] = (mean, sd) to compare against.
    Returns the per-gene details and the overall panel score.
    """
    genes = []
    for p in panel:
        gene = p["gene"]
        if gene not in sample["values"] or gene not in centers:
            continue
        value = sample["values"][gene]
        mean, sd = centers[gene]
        z = (value - mean) / sd
        contribution = p["sign"] * z
        genes.append({
            "gene": gene,
            "value": value,
            "z": z,
            "sign": p["sign"],
            "contribution": contribution,
        })

    if len(genes) == 0:
        return {"genes": [], "found": 0, "score": None, "mean_abs_z": 0}

    score = average([g["contribution"] for g in genes])
    mean_abs_z = average([abs(g["z"]) for g in genes])
    return {"genes": genes, "found": len(genes), "score": score, "mean_abs_z": mean_abs_z}


def group_label(text):
    """
    Turn the text in a 'group' column into 1 (Alzheimer's), 0 (control)
    or None (unknown / not given).
    """
    t = text.strip().lower()
    if t == "":
        return None
    for word in ("control", "normal", "healthy", "ctl", "ctrl", "nc"):
        if t == word or t.startswith(word):
            return 0
    for word in ("alzheimer", "ad", "case", "disease", "patient"):
        if t == word or t.startswith(word):
            return 1
    return None


def auc(scores, labels):
    """
    Area under the ROC curve: the chance that a random Alzheimer's sample
    scores higher than a random control sample. 0.5 = coin flip, 1 = perfect.
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


def analyze(panel, samples, cutoff, mode="reference"):
    """
    Score every sample and decide which pattern it looks like.
    mode is "reference" or "batch" (see the comment above).
    Returns a dictionary with the results and a few summary numbers.
    """
    needed = math.ceil(len(panel) * MIN_GENE_FRACTION)

    if mode == "batch":
        if len(samples) < MIN_BATCH_SIZE:
            raise AnalysisError(
                "Batch mode needs at least " + str(MIN_BATCH_SIZE) + " samples, but your file has "
                + str(len(samples)) + ". Upload more samples or use the GSE48350 reference instead."
            )
        centers = batch_centers(panel, samples)
    else:
        centers = reference_centers(panel)

    results = []
    skipped = []
    for sample in samples:
        result = score_sample(sample, panel, centers)
        if result["found"] < needed:
            skipped.append(sample["name"])
            continue

        if result["score"] > cutoff:
            pattern = "Alzheimer's-like"
            pattern_class = "ad"
        else:
            pattern = "Control-like"
            pattern_class = "ctl"

        flags = []
        if result["found"] < len(panel):
            missing = [p["gene"] for p in panel if p["gene"] not in sample["values"]]
            flags.append("missing " + ", ".join(missing))
        if mode == "reference" and result["mean_abs_z"] > 5:
            flags.append("scale looks off")

        biggest = max(abs(g["contribution"]) for g in result["genes"])

        results.append({
            "name": sample["name"],
            "group": sample["group"],
            "label": group_label(sample["group"]),
            "found": result["found"],
            "score": result["score"],
            "genes": result["genes"],
            "pattern": pattern,
            "pattern_class": pattern_class,
            "flags": flags,
            "biggest": biggest,
        })

    if len(results) == 0:
        raise AnalysisError(
            "No sample contains enough panel genes (each sample needs at least "
            + str(needed) + " of the " + str(len(panel)) + " genes)."
        )

    # Same bar scale for all samples so they can be compared
    bar_scale = max(3.0, max(r["biggest"] for r in results))
    for r in results:
        for g in r["genes"]:
            g["bar_width"] = min(100, abs(g["contribution"]) / bar_scale * 100)

    scale_warnings = len([r for r in results if "scale looks off" in r["flags"]])
    ad_count = len([r for r in results if r["pattern_class"] == "ad"])
    return {
        "results": results,
        "skipped": skipped,
        "needed": needed,
        "ad_count": ad_count,
        "ctl_count": len(results) - ad_count,
        "cutoff": cutoff,
        "mode": mode,
        "panel_size": len(panel),
        "scale_warnings": scale_warnings,
        "check": check_against_labels(results, panel, cutoff),
    }


def check_against_labels(results, panel, cutoff):
    """
    If the user told us which samples are Alzheimer's and which are controls,
    measure how well the panel separates them. Returns None if there are not
    at least 3 samples in each group.
    """
    labelled = [r for r in results if r["label"] is not None]
    labels = [r["label"] for r in labelled]
    n_ad = labels.count(1)
    n_ctl = labels.count(0)
    if n_ad < 3 or n_ctl < 3:
        return None

    scores = [r["score"] for r in labelled]
    true_pos = len([r for r in labelled if r["label"] == 1 and r["score"] > cutoff])
    true_neg = len([r for r in labelled if r["label"] == 0 and r["score"] <= cutoff])

    # Does each gene move in the direction the panel expects?
    genes = []
    for p in panel:
        ad_values = [r_value(r, p["gene"]) for r in labelled if r["label"] == 1]
        ctl_values = [r_value(r, p["gene"]) for r in labelled if r["label"] == 0]
        ad_values = [v for v in ad_values if v is not None]
        ctl_values = [v for v in ctl_values if v is not None]
        if len(ad_values) == 0 or len(ctl_values) == 0:
            continue
        difference = average(ad_values) - average(ctl_values)
        genes.append({
            "gene": p["gene"],
            "expected": "higher" if p["sign"] > 0 else "lower",
            "observed": "higher" if difference > 0 else "lower",
            "agrees": (difference > 0) == (p["sign"] > 0),
            "difference": difference,
        })

    return {
        "n_ad": n_ad,
        "n_ctl": n_ctl,
        "auc": auc(scores, labels),
        "accuracy": (true_pos + true_neg) / len(labelled),
        "correct": true_pos + true_neg,
        "total": len(labelled),
        "sensitivity": true_pos / n_ad,
        "specificity": true_neg / n_ctl,
        "genes": genes,
        "genes_agree": len([g for g in genes if g["agrees"]]),
        "roc": roc_svg(roc_curve(scores, labels)),
    }


def r_value(result, gene):
    """The measured value of one gene in one scored sample (or None)."""
    for g in result["genes"]:
        if g["gene"] == gene:
            return g["value"]
    return None


def describe_auc(value):
    """Plain-words meaning of an AUC value."""
    if value >= 0.9:
        return "Excellent separation"
    if value >= 0.8:
        return "Good separation"
    if value >= 0.7:
        return "Fair separation"
    if value > 0.6:
        return "Weak separation"
    if value >= 0.4:
        return "No real separation (close to a coin flip)"
    return "Reversed: Alzheimer's samples tend to score lower than controls"


# ---------------------------------------------------------------------------
# The "Where do your samples fall?" chart
# ---------------------------------------------------------------------------

def build_chart(report, reference=None):
    """
    Work out the x position (in pixels) of every dot on the chart.
    The HTML template then draws them as an SVG picture.
    """
    width = 720
    left = 30
    right = width - 30

    results = report["results"]
    all_values = [r["score"] for r in results] + [report["cutoff"]]
    if reference:
        all_values += reference["controls"] + reference["alzheimers"]

    low = min(all_values)
    high = max(all_values)
    if low == high:
        low -= 1
        high += 1
    padding = (high - low) * 0.1
    low -= padding
    high += padding

    def x_of(value):
        return left + (value - low) / (high - low) * (right - left)

    def points(items):
        return [{"x": x_of(r["score"]), "label": r["name"] + ": %.2f" % r["score"]} for r in items]

    # Rows of dots, from top to bottom
    rows = []
    if reference:
        rows.append({"label": "Reference controls (GSE48350)", "kind": "ctl",
                     "points": [{"x": x_of(v), "label": "%.2f" % v} for v in reference["controls"]]})
        rows.append({"label": "Reference Alzheimer's (GSE48350)", "kind": "ad",
                     "points": [{"x": x_of(v), "label": "%.2f" % v} for v in reference["alzheimers"]]})

    if report["check"]:
        # The user gave groups: draw their controls and Alzheimer's separately
        rows.append({"label": "Your samples labelled control", "kind": "ctl-you",
                     "points": points([r for r in results if r["label"] == 0])})
        rows.append({"label": "Your samples labelled Alzheimer's", "kind": "ad-you",
                     "points": points([r for r in results if r["label"] == 1])})
        unlabelled = [r for r in results if r["label"] is None]
        if unlabelled:
            rows.append({"label": "Your samples without a label", "kind": "you", "points": points(unlabelled)})
    else:
        rows.append({"label": "Your samples", "kind": "you", "points": points(results)})

    # Give every row its own height and spread dots a little up/down (jitter)
    row_height = 56
    top = 40
    for n, row in enumerate(rows):
        row["y"] = top + n * row_height
        for i, point in enumerate(row["points"]):
            point["y"] = row["y"] + ((i * 7) % 5 - 2) * 4

    bottom = top + len(rows) * row_height + 8   # room for the zone labels

    ticks = []
    for i in range(5):
        value = low + (high - low) * i / 4
        ticks.append({"x": x_of(value), "label": "%.1f" % value})

    return {
        "width": width,
        "height": bottom + 40,
        "left": left,
        "right": right,
        "top": top - 26,
        "axis_y": bottom,
        "rows": rows,
        "ticks": ticks,
        "cutoff_x": x_of(report["cutoff"]),
        "cutoff": report["cutoff"],
    }


# ---------------------------------------------------------------------------
# Results of the independent check on GSE5281 (made by validation/validate_gse5281.py)
# ---------------------------------------------------------------------------

VALIDATION_FOLDER = os.path.join(os.path.dirname(__file__), "validation", "results")
VALIDATION_SUMMARY = os.path.join(VALIDATION_FOLDER, "gse5281_summary.json")
VALIDATION_BATCH = os.path.join(VALIDATION_FOLDER, "gse5281_panel_genes.csv")


def load_validation():
    """Read the GSE5281 summary, or return None if the script was not run yet."""
    import json
    if not os.path.exists(VALIDATION_SUMMARY):
        return None
    with open(VALIDATION_SUMMARY) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Helpers for the pictures on the website
# ---------------------------------------------------------------------------

def roc_curve(scores, labels):
    """
    Points of the ROC curve as (false positive rate, true positive rate).
    We walk down the scores from highest to lowest; every step lowers the
    cutoff a little, so more samples get called Alzheimer's-like.
    """
    pairs = sorted(zip(scores, labels), key=lambda pair: -pair[0])
    positives = labels.count(1)
    negatives = labels.count(0)
    points = [(0.0, 0.0)]
    true_pos = 0
    false_pos = 0
    i = 0
    while i < len(pairs):
        current = pairs[i][0]
        # samples with exactly the same score move together
        while i < len(pairs) and pairs[i][0] == current:
            if pairs[i][1] == 1:
                true_pos += 1
            else:
                false_pos += 1
            i += 1
        points.append((false_pos / negatives, true_pos / positives))
    return points


def roc_svg(points, size=220, pad=34):
    """Turn ROC points into an SVG 'points' string for a <polyline>."""
    coords = []
    for fpr, tpr in points:
        x = pad + fpr * size
        y = pad + size - tpr * size
        coords.append("%.1f,%.1f" % (x, y))
    return " ".join(coords)


def heat_style(value, limit=2.5):
    """
    Colour for one heatmap cell. Positive = toward Alzheimer's (coral),
    negative = toward control (teal). Stronger values = stronger colour.
    Returns a CSS style string.
    """
    strength = min(1.0, abs(value) / limit)
    percent = round(strength * 85)
    colour = "--ad" if value > 0 else "--ctl"
    text = "#fff" if percent > 50 else "var(--ink)"
    return "background: color-mix(in srgb, var(%s) %d%%, var(--cell)); color: %s;" % (colour, percent, text)


def validation_roc():
    """ROC curve of the GSE5281 validation, read from its scores file."""
    path = os.path.join(VALIDATION_FOLDER, "gse5281_scores.csv")
    if not os.path.exists(path):
        return None
    scores = []
    labels = []
    with open(path) as f:
        for row in csv.DictReader(f):
            scores.append(float(row["panel_score"]))
            labels.append(group_label(row["group"]))
    return roc_svg(roc_curve(scores, labels))
