"""
app.py
------
The website part of the project, built with Flask.

Run it with:
    python app.py
and open http://127.0.0.1:5000 in your browser.

Pages:
    /                  -> home page with the upload form
    /analyze           -> receives the form, runs the analysis, shows the results
    /example           -> made-up example data (reference mode)
    /example/gse5281   -> real data: 23 hippocampus samples from GSE5281 (batch mode)
    /validation        -> how the panel did on the independent GSE5281 dataset
    /about             -> how the panel score works
"""

from flask import Flask, render_template, request

import analysis

app = Flask(__name__)

# Whole-genome expression files can be big, so allow uploads up to 50 MB.
app.config["MAX_CONTENT_LENGTH"] = 50 * 1024 * 1024

REFERENCE_CUTOFF_TEXT = "95th percentile of the 18 GSE48350 control scores (measured on the same data the panel was chosen from)"
BATCH_CUTOFF_TEXT = "0 = the average of your batch"


def get_text(file_field, text_field):
    """
    The user can either upload a file OR paste text.
    If a file was uploaded we use it, otherwise we use the pasted text.
    """
    uploaded = request.files.get(file_field)
    if uploaded and uploaded.filename:
        # "utf-8-sig" also removes the invisible marker Excel adds to CSV files
        return uploaded.read().decode("utf-8-sig", errors="replace")
    return request.form.get(text_field, "").strip()


def show_form(error=None, form=None):
    """Draw the home page. Used on first visit and when something went wrong."""
    form = form or {}
    panel_text = form.get("panel_text") or analysis.default_panel_text()
    return render_template(
        "index.html",
        error=error,
        form=form,
        panel_text=panel_text,
        panel=analysis.parse_panel(analysis.default_panel_text()),
        default_cutoff=analysis.DEFAULT_CUTOFF,
        min_batch=analysis.MIN_BATCH_SIZE,
        validation=analysis.load_validation(),
    )


def show_results(panel, samples, cutoff, cutoff_source, mode, reference=None, title=None):
    """Score the samples and draw the results page."""
    report = analysis.analyze(panel, samples, cutoff, mode)
    chart = analysis.build_chart(report, reference)
    auc_words = analysis.describe_auc(report["check"]["auc"]) if report["check"] else None
    return render_template(
        "results.html",
        report=report,
        chart=chart,
        cutoff_source=cutoff_source,
        reference=reference,
        auc_words=auc_words,
        title=title,
        min_batch=analysis.MIN_BATCH_SIZE,
    )


@app.route("/")
def index():
    return show_form()


@app.route("/analyze", methods=["POST"])
def analyze():
    try:
        # 1. Panel reference (use the built-in one if the box is empty)
        panel_text = get_text("panel_file", "panel_text")
        if panel_text == "":
            panel_text = analysis.default_panel_text()
        panel = analysis.parse_panel(panel_text)

        # 2. Expression data
        expression_text = get_text("expression_file", "expression_text")
        if expression_text == "":
            raise analysis.AnalysisError("Please upload or paste your gene expression data first.")
        samples = analysis.parse_expression(expression_text, panel)

        # 3. What to compare against: the GSE48350 numbers, or the batch itself
        mode = request.form.get("compare_mode", "auto")
        if mode == "auto":
            mode = "batch" if len(samples) >= analysis.MIN_BATCH_SIZE else "reference"

        # 4. Cutoff
        cutoff_mode = request.form.get("cutoff_mode", "default")
        reference = None
        if cutoff_mode == "custom":
            typed = request.form.get("cutoff_value", "").strip()
            if not analysis.is_number(typed):
                raise analysis.AnalysisError("Please type a number for the custom cutoff.")
            cutoff = float(typed)
            cutoff_source = "custom value you typed"
        elif cutoff_mode == "reference" and mode == "reference":
            scores_text = get_text("scores_file", "scores_text")
            if scores_text == "":
                raise analysis.AnalysisError(
                    "You chose 'from reference scores' but did not add a scores file."
                )
            reference = analysis.parse_reference_scores(scores_text)
            cutoff = reference["cutoff"]
            cutoff_source = "95th percentile of your reference control scores"
        elif mode == "batch":
            cutoff = 0.0
            cutoff_source = BATCH_CUTOFF_TEXT
        else:
            cutoff = analysis.DEFAULT_CUTOFF
            cutoff_source = REFERENCE_CUTOFF_TEXT

        return show_results(panel, samples, cutoff, cutoff_source, mode, reference)

    except analysis.AnalysisError as error:
        # Show the form again with the error and keep what the user typed
        return show_form(error=str(error), form=request.form)


@app.route("/example")
def example():
    panel = analysis.parse_panel(analysis.default_panel_text())
    samples = analysis.parse_expression(analysis.read_file_text(analysis.EXAMPLE_FILE), panel)
    return show_results(panel, samples, analysis.DEFAULT_CUTOFF, REFERENCE_CUTOFF_TEXT,
                        "reference", title="Made-up example (5 samples)")


@app.route("/example/gse5281")
def example_gse5281():
    panel = analysis.parse_panel(analysis.default_panel_text())
    samples = analysis.parse_expression(analysis.read_file_text(analysis.VALIDATION_BATCH), panel)
    return show_results(panel, samples, 0.0, BATCH_CUTOFF_TEXT, "batch",
                        title="GSE5281 hippocampus (23 real samples)")


@app.route("/validation")
def validation():
    return render_template("validation.html", v=analysis.load_validation())


@app.route("/about")
def about():
    panel = analysis.parse_panel(analysis.default_panel_text())
    return render_template("about.html", panel=panel, default_cutoff=analysis.DEFAULT_CUTOFF,
                           min_batch=analysis.MIN_BATCH_SIZE)


@app.errorhandler(413)
def file_too_large(error):
    return show_form(error="That file is larger than 50 MB. Try keeping only the gene rows you need.")


if __name__ == "__main__":
    app.run(debug=True)
