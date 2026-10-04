"""Tests for Module C: the mock DE generator and the volcano-plot CLI."""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def de_results() -> pd.DataFrame:
    """DESeq2-style results (note DESeq2's lower-case 'log2FoldChange')."""
    return pd.DataFrame(
        {
            "gene": ["UP1", "DOWN1", "SMALL_FC", "NOT_SIG", "NO_P", "ZERO_P"],
            "log2FoldChange": [2.0, -1.5, 0.5, 3.0, 1.0, 4.0],
            "pvalue": [1e-3, 1e-2, 1e-10, 0.2, np.nan, 0.0],
        }
    )


def test_column_matching_ignores_case(plot_volcano, de_results):
    assert plot_volcano.find_column(de_results, "Log2FoldChange") == "log2FoldChange"
    with pytest.raises(ValueError, match="not found"):
        plot_volcano.find_column(de_results, "logFC")


def test_gene_column_auto_detection(plot_volcano):
    r_export = pd.DataFrame({"Unnamed: 0": ["TP53"], "log2FoldChange": [1.0]})
    numeric_only = pd.DataFrame({"log2FoldChange": [1.0], "pvalue": [0.1]})
    assert plot_volcano.find_gene_column(r_export, None) == "Unnamed: 0"
    assert plot_volcano.find_gene_column(numeric_only, None) is None


def test_prepare_data_classifies_and_cleans(plot_volcano, de_results):
    data = plot_volcano.prepare_data(
        de_results, "log2FoldChange", "pvalue", "gene", 1.0, 0.05
    )
    regulation = dict(zip(data["gene"], data["regulation"]))

    assert "NO_P" not in regulation  # missing p-value is skipped
    assert regulation["UP1"] == plot_volcano.UP
    assert regulation["DOWN1"] == plot_volcano.DOWN
    assert regulation["SMALL_FC"] == plot_volcano.NOT_SIGNIFICANT
    assert regulation["NOT_SIG"] == plot_volcano.NOT_SIGNIFICANT
    assert np.isfinite(data["neg_log10_p"]).all()  # p = 0 is clipped, not infinite


def test_prepare_data_rejects_values_that_are_not_p_values(plot_volcano, de_results):
    with pytest.raises(ValueError, match="outside 0-1"):
        plot_volcano.prepare_data(
            de_results, "log2FoldChange", "log2FoldChange", "gene", 1.0, 0.05
        )


def test_labels_are_split_between_directions(plot_volcano):
    data = pd.DataFrame(
        {
            "gene": ["U1", "U2", "U3", "D1", "D2"],
            "pvalue": [1e-5, 1e-4, 1e-3, 1e-6, 1e-2],
            "regulation": [plot_volcano.UP] * 3 + [plot_volcano.DOWN] * 2,
        }
    )
    labels = plot_volcano.pick_labels(data, top_n=4)
    assert sorted(labels["gene"]) == ["D1", "D2", "U1", "U2"]


def test_cli_creates_png(plot_volcano, de_results, tmp_path):
    source = tmp_path / "deseq2_results.csv"
    de_results.to_csv(source, index=False)
    image = tmp_path / "plots" / "volcano.png"

    exit_code = plot_volcano.main(["--input", str(source), "--output", str(image)])

    assert exit_code == 0
    assert image.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_cli_reports_missing_column(plot_volcano, de_results, tmp_path):
    source = tmp_path / "results.csv"
    de_results.drop(columns="pvalue").to_csv(source, index=False)
    exit_code = plot_volcano.main(
        ["--input", str(source), "--output", str(tmp_path / "volcano.png")]
    )
    assert exit_code == 1


def test_benjamini_hochberg(generate_mock_de_results):
    pvalues = np.array([0.01, 0.04, 0.03, 0.005, np.nan])
    adjusted = generate_mock_de_results.benjamini_hochberg(pvalues)
    np.testing.assert_allclose(adjusted[:4], [0.02, 0.04, 0.04, 0.02])
    assert np.isnan(adjusted[4])


def test_mock_generator_writes_deseq2_style_table(generate_mock_de_results, tmp_path):
    output = tmp_path / "mock.csv"
    assert generate_mock_de_results.main(["-o", str(output), "--n-genes", "600"]) == 0

    table = pd.read_csv(output)
    assert list(table.columns) == [
        "gene",
        "baseMean",
        "Log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
    ]
    assert len(table) == 600
    both = table.dropna(subset=["pvalue", "padj"])
    assert (both["padj"] >= both["pvalue"]).all()  # FDR never shrinks a p-value
