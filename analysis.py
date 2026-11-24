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
  3. For every gene, compute a z-score:  z = (value - control_mean) / control_sd
     and multiply it by the gene's sign (+1 = higher in Alzheimer's,
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
        samples.append({"name": names[n], "group": "", "values": values})
    return samples


# ---------------------------------------------------------------------------
# Step 3 and 4: scoring
# ---------------------------------------------------------------------------

def score_sample(sample, panel):
    """
    Compare one sample with the panel.
    Returns the per-gene details and the overall panel score.
    """
    genes = []
    for p in panel:
        gene = p["gene"]
        if gene not in sample["values"]:
            continue
        value = sample["values"][gene]
        z = (value - p["mean"]) / p["sd"]
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


def analyze(panel, samples, cutoff):
    """
    Score every sample and decide which pattern it looks like.
    Returns a dictionary with the results and a few summary numbers.
    """
    needed = math.ceil(len(panel) * MIN_GENE_FRACTION)
    results = []
    skipped = []

    for sample in samples:
        result = score_sample(sample, panel)
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
        if result["mean_abs_z"] > 5:
            flags.append("scale looks off")

        # Largest absolute contribution, used to draw the bars at the same scale
        biggest = max(abs(g["contribution"]) for g in result["genes"])

        results.append({
            "name": sample["name"],
            "group": sample["group"],
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

    ad_count = len([r for r in results if r["pattern_class"] == "ad"])
    return {
        "results": results,
        "skipped": skipped,
        "needed": needed,
        "ad_count": ad_count,
        "ctl_count": len(results) - ad_count,
        "cutoff": cutoff,
        "panel_size": len(panel),
    }


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

    sample_scores = [r["score"] for r in report["results"]]
    all_values = sample_scores + [report["cutoff"]]
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

    # Rows of dots, from top to bottom
    rows = []
    if reference:
        rows.append({
            "label": "Reference controls",
            "kind": "ctl",
            "points": [{"x": x_of(v), "label": "%.2f" % v} for v in reference["controls"]],
        })
        rows.append({
            "label": "Reference Alzheimer's",
            "kind": "ad",
            "points": [{"x": x_of(v), "label": "%.2f" % v} for v in reference["alzheimers"]],
        })
    rows.append({
        "label": "Your samples",
        "kind": "you",
        "points": [
            {"x": x_of(r["score"]), "label": r["name"] + ": %.2f" % r["score"]}
            for r in report["results"]
        ],
    })

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
