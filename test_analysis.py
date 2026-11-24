"""
Simple checks for analysis.py. Run with:   python test_analysis.py
Each test uses plain `assert`: if something is wrong, Python stops with an error.
"""

import analysis


def test_percentile():
    assert analysis.percentile([1, 2, 3, 4, 5], 0.5) == 3
    assert analysis.percentile([0, 10], 0.95) == 9.5


def test_panel_has_seven_genes():
    panel = analysis.parse_panel(analysis.default_panel_text())
    assert len(panel) == 7
    assert panel[0]["gene"] == "ANKIB1"
    assert panel[1]["sign"] == -1


def test_sample_equal_to_control_mean_scores_zero():
    panel = analysis.parse_panel(analysis.default_panel_text())
    values = {p["gene"]: p["mean"] for p in panel}
    centers = analysis.reference_centers(panel)
    result = analysis.score_sample({"name": "x", "group": "", "values": values}, panel, centers)
    assert abs(result["score"]) < 1e-9


def test_both_layouts_give_the_same_score():
    panel = analysis.parse_panel(analysis.default_panel_text())
    genes_as_rows = "gene,S1\n" + "\n".join(p["gene"] + "," + str(p["mean"] + p["sign"] * p["sd"]) for p in panel)
    samples_as_rows = ("sample," + ",".join(p["gene"] for p in panel) + "\nS1,"
                       + ",".join(str(p["mean"] + p["sign"] * p["sd"]) for p in panel))
    a = analysis.parse_expression(genes_as_rows, panel)
    b = analysis.parse_expression(samples_as_rows, panel)
    centers = analysis.reference_centers(panel)
    score_a = analysis.score_sample(a[0], panel, centers)["score"]
    score_b = analysis.score_sample(b[0], panel, centers)["score"]
    assert abs(score_a - 1.0) < 1e-9   # every gene is 1 SD toward Alzheimer's
    assert abs(score_a - score_b) < 1e-9


def test_example_file():
    panel = analysis.parse_panel(analysis.default_panel_text())
    samples = analysis.parse_expression(analysis.read_file_text(analysis.EXAMPLE_FILE), panel)
    report = analysis.analyze(panel, samples, analysis.DEFAULT_CUTOFF)
    assert len(report["results"]) == 5
    assert report["ad_count"] == 2


def test_no_panel_genes_gives_friendly_error():
    panel = analysis.parse_panel(analysis.default_panel_text())
    try:
        analysis.parse_expression("gene,S1\nAPP,1.0\nMAPT,2.0", panel)
    except analysis.AnalysisError:
        return
    assert False, "expected an AnalysisError"


def test_group_labels():
    assert analysis.group_label("control") == 0
    assert analysis.group_label("Normal") == 0
    assert analysis.group_label("AD") == 1
    assert analysis.group_label("Alzheimer's Disease") == 1
    assert analysis.group_label("") is None


def test_auc():
    assert analysis.auc([1, 2, 3, 4], [0, 0, 1, 1]) == 1.0   # perfect
    assert analysis.auc([4, 3, 2, 1], [0, 0, 1, 1]) == 0.0   # reversed
    assert analysis.auc([1, 1, 1, 1], [0, 0, 1, 1]) == 0.5   # all tied


def test_batch_mode_needs_enough_samples():
    panel = analysis.parse_panel(analysis.default_panel_text())
    samples = analysis.parse_expression(analysis.read_file_text(analysis.EXAMPLE_FILE), panel)
    try:
        analysis.analyze(panel, samples, 0.0, "batch")
    except analysis.AnalysisError:
        return
    assert False, "expected an AnalysisError for a batch of 5"


def test_batch_mode_ignores_scale():
    # Adding 10 to every value (a different lab's scale) must not change batch scores
    panel = analysis.parse_panel(analysis.default_panel_text())
    text = analysis.read_file_text(analysis.VALIDATION_BATCH)
    samples = analysis.parse_expression(text, panel)
    shifted = analysis.parse_expression(text, panel)
    for s in shifted:
        for g in s["values"]:
            s["values"][g] += 10
    a = analysis.analyze(panel, samples, 0.0, "batch")
    b = analysis.analyze(panel, shifted, 0.0, "batch")
    for r1, r2 in zip(a["results"], b["results"]):
        assert abs(r1["score"] - r2["score"]) < 1e-9


def test_gse5281_matches_validation_script():
    # The website and validation/validate_gse5281.py must give the same AUC
    panel = analysis.parse_panel(analysis.default_panel_text())
    samples = analysis.parse_expression(analysis.read_file_text(analysis.VALIDATION_BATCH), panel)
    report = analysis.analyze(panel, samples, 0.0, "batch")
    assert report["check"]["n_ctl"] == 13 and report["check"]["n_ad"] == 10
    assert abs(report["check"]["auc"] - analysis.load_validation()["auc"]) < 0.001


def test_group_row_in_genes_as_rows_layout():
    panel = analysis.parse_panel(analysis.default_panel_text())
    text = "gene,S1,S2\ngroup,control,AD\n" + "\n".join(p["gene"] + ",1.0,2.0" for p in panel)
    samples = analysis.parse_expression(text, panel)
    assert [s["group"] for s in samples] == ["control", "AD"]
    assert len(samples[0]["values"]) == 7


if __name__ == "__main__":
    tests = [value for name, value in list(globals().items()) if name.startswith("test_")]
    for test in tests:
        test()
        print("OK  ", test.__name__)
    print("All", len(tests), "tests passed.")
