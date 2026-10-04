"""Tests for Module A: the messy-data generator and the cleaning CLI."""

import pandas as pd
import pytest


def make_raw() -> pd.DataFrame:
    """A tiny spreadsheet with one example of every problem the cleaner fixes.

    Rows: TP53 | an exact duplicate once normalised | BRCA1 twice with different
    counts | MYC with a negative count | no symbol | EGFR with a text
    placeholder | a completely empty row. One header has a trailing space.
    """
    return pd.DataFrame(
        {
            "Gene_Symbol": [
                "TP53",
                " tp53",
                "brca1",
                "BRCA1 ",
                "MYC",
                None,
                "EGFR",
                None,
            ],
            "Control_1": [10, 10, 5, 7, 8, 3, "n.d.", None],
            "Treated_1 ": [20, 20, 6, 9, -4, 3, 11, None],
        }
    )


def test_normalize_symbols(clean_data):
    symbols = pd.Series([" tp53", "Brca1 ", "", "  ", None])
    result = clean_data.normalize_symbols(symbols)
    assert result.iloc[:2].tolist() == ["TP53", "BRCA1"]
    assert result.iloc[2:].isna().all()


def test_default_pipeline_fixes_every_problem(clean_data):
    clean, report = clean_data.clean_counts(make_raw())

    expected = pd.DataFrame(
        {"Gene_Symbol": ["BRCA1", "TP53"], "Control_1": [12, 10], "Treated_1": [15, 20]}
    )
    pd.testing.assert_frame_equal(clean, expected, check_dtype=False)
    assert (report.empty_rows, report.missing_symbol_rows) == (1, 1)
    assert report.invalid_cells == 2  # "n.d." and the negative count
    assert report.exact_duplicates == 1
    assert report.rows_with_missing == 2  # MYC and EGFR are dropped
    assert report.merged_duplicates == 1  # BRCA1 summed
    assert report.output_genes == 2


@pytest.mark.parametrize(
    ("strategy", "myc", "egfr"),
    [("zero", [8, 0], [0, 11]), ("median", [8, 8], [11, 11])],
)
def test_missing_value_strategies_keep_genes(clean_data, strategy, myc, egfr):
    clean, _ = clean_data.clean_counts(make_raw(), na_strategy=strategy)
    counts = clean.set_index("Gene_Symbol")
    assert counts.loc["MYC"].tolist() == myc
    assert counts.loc["EGFR"].tolist() == egfr


@pytest.mark.parametrize(("strategy", "brca1"), [("mean", [6, 8]), ("first", [5, 6])])
def test_duplicate_strategies(clean_data, strategy, brca1):
    clean, _ = clean_data.clean_counts(make_raw(), duplicate_strategy=strategy)
    assert clean.set_index("Gene_Symbol").loc["BRCA1"].tolist() == brca1


def test_unknown_gene_column_raises(clean_data):
    with pytest.raises(ValueError, match="not found"):
        clean_data.clean_counts(make_raw(), gene_col="Gene_ID")


def test_cli_end_to_end(clean_data, tmp_path):
    source = tmp_path / "raw.xlsx"
    make_raw().to_excel(source, index=False)
    target = tmp_path / "out" / "clean.csv"

    exit_code = clean_data.main(["--input", str(source), "--output", str(target)])

    assert exit_code == 0
    assert pd.read_csv(target)["Gene_Symbol"].tolist() == ["BRCA1", "TP53"]


def test_cli_reports_missing_input(clean_data, tmp_path):
    exit_code = clean_data.main(
        ["--input", str(tmp_path / "nope.xlsx"), "--output", str(tmp_path / "x.csv")]
    )
    assert exit_code == 1


def test_cli_refuses_to_overwrite_input(clean_data, tmp_path):
    source = tmp_path / "raw.csv"
    with pytest.raises(SystemExit) as excinfo:
        clean_data.parse_args(["--input", str(source), "--output", str(source)])
    assert excinfo.value.code == 2


def test_generated_workbook_is_messy_but_cleanable(
    generate_mock_rnaseq, clean_data, tmp_path
):
    workbook = tmp_path / "raw_counts.xlsx"
    assert generate_mock_rnaseq.main(["--output", str(workbook), "--seed", "1"]) == 0

    raw = pd.read_excel(workbook)
    assert raw["Gene_Symbol"].duplicated().any()
    assert raw.isna().to_numpy().any()

    clean, report = clean_data.clean_counts(raw)
    assert report.exact_duplicates > 0
    assert not clean.isna().to_numpy().any()
    assert clean["Gene_Symbol"].is_unique
    assert (clean["Gene_Symbol"] == clean["Gene_Symbol"].str.upper()).all()
