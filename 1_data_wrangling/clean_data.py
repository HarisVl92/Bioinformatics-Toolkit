#!/usr/bin/env python3
"""Clean a messy RNA-seq raw-count table into an analysis-ready count matrix.

Fixes the problems typically found in hand-edited spreadsheets:

* gene symbols in mixed case or padded with whitespace
* empty rows, rows without a gene symbol and exact duplicate rows
* text placeholders ("NA", "n/a", "-") and negative numbers in count columns
* missing counts (configurable: drop the gene, fill with 0, or fill with the
  gene's median across samples)
* gene symbols reported more than once (configurable: sum, mean or first)

Every step is logged and summarised, so the cleaning is fully auditable.
The input file is never modified.

Example:
    python clean_data.py --input data/raw_counts.xlsx --output output/clean_counts.csv
"""

import argparse
import logging
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

LOGGER = logging.getLogger("clean_data")

EXCEL_SUFFIXES = {".xlsx", ".xlsm"}
TEXT_SEPARATORS = {".csv": ",", ".tsv": "\t", ".txt": "\t"}
NA_STRATEGIES = ("drop", "zero", "median")
DUPLICATE_STRATEGIES = ("sum", "mean", "first")


@dataclass
class CleaningReport:
    """What each cleaning step changed, for an auditable summary."""

    na_strategy: str
    duplicate_strategy: str
    input_rows: int = 0
    empty_rows: int = 0
    missing_symbol_rows: int = 0
    invalid_cells: int = 0
    exact_duplicates: int = 0
    rows_with_missing: int = 0
    merged_duplicates: int = 0
    output_genes: int = 0

    def log(self) -> None:
        """Write the summary table to the log."""
        na_action = {
            "drop": "removed",
            "zero": "filled with 0",
            "median": "filled with median",
        }[self.na_strategy]
        lines = [
            ("Input rows", self.input_rows),
            ("Empty rows removed", self.empty_rows),
            ("Rows without a gene symbol removed", self.missing_symbol_rows),
            ("Non-numeric/negative cells set to missing", self.invalid_cells),
            ("Exact duplicate rows removed", self.exact_duplicates),
            (f"Rows with missing counts {na_action}", self.rows_with_missing),
            (
                f"Duplicate symbols merged ({self.duplicate_strategy})",
                self.merged_duplicates,
            ),
            ("Genes in output", self.output_genes),
        ]
        LOGGER.info("Cleaning summary:")
        for label, value in lines:
            LOGGER.info("  %s %s %5d", label, "." * (44 - len(label)), value)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse and validate command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Clean a messy RNA-seq raw-count table (Excel, CSV or TSV)\n"
        "into a tidy, integer count matrix.",
        epilog="examples:\n"
        "  python clean_data.py --input data/raw_counts.xlsx "
        "--output output/clean_counts.csv\n"
        "  python clean_data.py -i counts.csv -o clean.tsv "
        "--na-strategy median --duplicate-strategy mean",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="raw count table (.xlsx, .xlsm, .csv or .tsv)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="destination file (.csv, or .tsv for tab-separated output)",
    )
    parser.add_argument(
        "--gene-col",
        default="Gene_Symbol",
        help="column holding gene symbols, case-insensitive (default: %(default)s)",
    )
    parser.add_argument(
        "--sheet",
        default="0",
        help="Excel sheet name or zero-based index (default: %(default)s)",
    )
    parser.add_argument(
        "--na-strategy",
        choices=NA_STRATEGIES,
        default="drop",
        help="how to handle missing counts: drop the gene, fill with 0, or fill "
        "with the gene's median across samples (default: %(default)s)",
    )
    parser.add_argument(
        "--duplicate-strategy",
        choices=DUPLICATE_STRATEGIES,
        default="sum",
        help="how to merge rows sharing a gene symbol (default: %(default)s)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="show debug messages"
    )
    args = parser.parse_args(argv)
    if args.output.resolve() == args.input.resolve():
        parser.error(
            "--output must differ from --input (raw data is never overwritten)"
        )
    return args


def configure_logging(verbose: bool) -> None:
    """Log to stderr as 'time | level | message'.

    Only this tool logs at INFO/DEBUG; third-party libraries report warnings only.
    """
    logging.basicConfig(
        level=logging.WARNING,
        format="%(asctime)s | %(levelname)-8s | %(message)s",
        datefmt="%H:%M:%S",
    )
    LOGGER.setLevel(logging.DEBUG if verbose else logging.INFO)


def load_table(path: Path, sheet: str | int = 0) -> pd.DataFrame:
    """Read an Excel workbook or a delimited text file."""
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path}")
    suffix = path.suffix.lower()
    if suffix in EXCEL_SUFFIXES:
        return pd.read_excel(path, sheet_name=sheet, engine="openpyxl")
    if suffix in TEXT_SEPARATORS:
        return pd.read_csv(path, sep=TEXT_SEPARATORS[suffix])
    raise ValueError(
        f"unsupported file type '{suffix}'; use .xlsx, .xlsm, .csv or .tsv"
    )


def find_column(columns: pd.Index, name: str) -> str:
    """Return the column called ``name``, ignoring case and surrounding spaces."""
    wanted = name.strip().lower()
    for column in columns:
        if str(column).strip().lower() == wanted:
            return column
    available = ", ".join(map(str, columns))
    raise ValueError(f"column '{name}' not found; available columns: {available}")


def select_sample_columns(frame: pd.DataFrame, gene_col: str) -> list[str]:
    """Return the count columns (mostly numeric); warn about the others."""
    samples, ignored = [], []
    for column in frame.columns.drop(gene_col):
        values = frame[column].dropna()
        numeric_share = pd.to_numeric(values, errors="coerce").notna().mean()
        if len(values) and numeric_share >= 0.5:
            samples.append(column)
        else:
            ignored.append(column)
    if ignored:
        LOGGER.warning("Ignoring non-numeric column(s): %s", ", ".join(ignored))
    if not samples:
        raise ValueError("no numeric sample columns found next to the gene column")
    return samples


def normalize_symbols(symbols: pd.Series) -> pd.Series:
    """Strip whitespace and upper-case gene symbols; blank symbols become missing.

    Upper case follows the HGNC convention for human genes (e.g. tp53 -> TP53).
    """
    cleaned = symbols.astype("string").str.strip().str.upper()
    return cleaned.replace("", pd.NA)


def coerce_counts(counts: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Convert count columns to numbers; text and negative values become NaN."""
    numeric = counts.apply(pd.to_numeric, errors="coerce")
    non_numeric = numeric.isna() & counts.notna()
    negative = numeric < 0

    if non_numeric.to_numpy().any():
        examples = sorted({str(v) for v in counts.to_numpy()[non_numeric.to_numpy()]})
        LOGGER.warning(
            "Set %d non-numeric cell(s) to missing: %s",
            non_numeric.to_numpy().sum(),
            ", ".join(repr(example) for example in examples[:5]),
        )
    if negative.to_numpy().any():
        LOGGER.warning("Set %d negative count(s) to missing", negative.to_numpy().sum())
    numeric = numeric.mask(negative)

    fractional = numeric.notna() & (numeric % 1 != 0)
    if fractional.to_numpy().any():
        LOGGER.warning(
            "%d value(s) are not whole numbers and will be rounded; raw counts "
            "are expected - is this normalised data (e.g. TPM)?",
            fractional.to_numpy().sum(),
        )
    invalid = int(non_numeric.to_numpy().sum() + negative.to_numpy().sum())
    return numeric, invalid


def handle_missing(
    frame: pd.DataFrame, sample_cols: list[str], strategy: str
) -> tuple[pd.DataFrame, int]:
    """Drop or impute genes that have missing counts."""
    has_missing = frame[sample_cols].isna().any(axis="columns")
    affected = int(has_missing.sum())
    if strategy == "drop":
        return frame.loc[~has_missing], affected
    if strategy == "zero":
        return frame.fillna({column: 0 for column in sample_cols}), affected

    # "median": the gene's median across the samples where it was measured.
    row_median = frame[sample_cols].median(axis="columns")
    filled = frame[sample_cols].apply(lambda column: column.fillna(row_median))
    frame = frame.drop(columns=sample_cols).join(filled)
    unobserved = frame[sample_cols].isna().any(axis="columns")
    if unobserved.any():
        LOGGER.warning(
            "Removed %d gene(s) with no observed counts at all", unobserved.sum()
        )
    return frame.loc[~unobserved], affected


def merge_duplicate_symbols(
    frame: pd.DataFrame, gene_col: str, sample_cols: list[str], strategy: str
) -> tuple[pd.DataFrame, int]:
    """Collapse rows that share a gene symbol into one row (sum, mean or first)."""
    duplicated = frame[gene_col].duplicated()
    if not duplicated.any():
        return frame, 0
    genes = sorted(frame.loc[duplicated, gene_col].unique())
    LOGGER.info("Merging duplicate gene symbols (%s): %s", strategy, ", ".join(genes))
    merged = frame.groupby(gene_col, sort=False, as_index=False)[sample_cols].agg(
        strategy
    )
    return merged, int(duplicated.sum())


def clean_counts(
    raw: pd.DataFrame,
    gene_col: str = "Gene_Symbol",
    na_strategy: str = "drop",
    duplicate_strategy: str = "sum",
) -> tuple[pd.DataFrame, CleaningReport]:
    """Run the full cleaning pipeline; return the clean matrix and a report."""
    report = CleaningReport(na_strategy, duplicate_strategy, input_rows=len(raw))

    frame = raw.rename(columns=lambda name: str(name).strip())
    empty_columns = [c for c in frame.columns if frame[c].isna().all()]
    if empty_columns:
        LOGGER.warning("Dropping empty column(s): %s", ", ".join(empty_columns))
        frame = frame.drop(columns=empty_columns)
    gene_col = find_column(frame.columns, gene_col)
    sample_cols = select_sample_columns(frame, gene_col)
    LOGGER.info("Sample columns (%d): %s", len(sample_cols), ", ".join(sample_cols))

    # 1. Gene symbols: normalise, then drop empty rows and rows without a symbol.
    frame = frame.assign(**{gene_col: normalize_symbols(frame[gene_col])})
    empty = frame[[gene_col, *sample_cols]].isna().all(axis="columns")
    report.empty_rows = int(empty.sum())
    frame = frame.loc[~empty]
    no_symbol = frame[gene_col].isna()
    report.missing_symbol_rows = int(no_symbol.sum())
    frame = frame.loc[~no_symbol]

    # 2. Counts: numbers only; text placeholders and negatives become missing.
    counts, report.invalid_cells = coerce_counts(frame[sample_cols])
    frame = pd.concat([frame[[gene_col]], counts], axis="columns")

    # 3. Exact duplicate rows (compared after normalisation, so " tp53" == "TP53").
    duplicated = frame.duplicated()
    report.exact_duplicates = int(duplicated.sum())
    frame = frame.loc[~duplicated]

    # 4. Missing counts, then 5. genes reported more than once.
    frame, report.rows_with_missing = handle_missing(frame, sample_cols, na_strategy)
    frame, report.merged_duplicates = merge_duplicate_symbols(
        frame, gene_col, sample_cols, duplicate_strategy
    )

    # 6. Integer counts, sorted by gene symbol.
    counts = frame[sample_cols].round().astype("int64")
    frame = pd.concat([frame[[gene_col]], counts], axis="columns")
    frame = frame.sort_values(gene_col, ignore_index=True)
    report.output_genes = len(frame)
    return frame, report


def write_table(frame: pd.DataFrame, path: Path) -> None:
    """Write CSV (or TSV for .tsv/.txt file names), creating folders as needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    separator = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    frame.to_csv(path, sep=separator, index=False)


def main(argv: list[str] | None = None) -> int:
    """Run the cleaning pipeline and return a process exit code."""
    args = parse_args(argv)
    configure_logging(args.verbose)
    sheet = int(args.sheet) if args.sheet.isdigit() else args.sheet

    try:
        raw = load_table(args.input, sheet)
        LOGGER.info("Loaded %d rows x %d columns from %s", *raw.shape, args.input)
        clean, report = clean_counts(
            raw, args.gene_col, args.na_strategy, args.duplicate_strategy
        )
        write_table(clean, args.output)
    except FileNotFoundError as exc:
        LOGGER.error("%s", exc)
        return 1
    except (ValueError, IndexError) as exc:
        LOGGER.error("Invalid input: %s", exc)
        return 1
    except zipfile.BadZipFile:
        LOGGER.error("%s is not a valid Excel workbook (corrupted file?)", args.input)
        return 1
    except OSError as exc:
        LOGGER.error("File error: %s", exc)
        return 1
    except KeyboardInterrupt:
        LOGGER.error("Interrupted by user")
        return 130

    report.log()
    LOGGER.info(
        "Saved %d genes x %d samples to %s",
        len(clean),
        clean.shape[1] - 1,
        args.output,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
