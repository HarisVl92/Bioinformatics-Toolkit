#!/usr/bin/env python3
"""Generate a realistic, deliberately messy RNA-seq raw-count spreadsheet.

The workbook mimics a count matrix that was exported to Excel and then edited
by hand. It contains the problems that ``clean_data.py`` is built to fix:

* gene symbols in mixed case or padded with spaces (``tp53``, `` Brca1``)
* exact duplicate rows, re-formatted copies of rows, and gene symbols that
  appear twice with different counts (e.g. two features of the same gene)
* blank cells and text placeholders (``NA``, ``n/a``, ``-``, ``n.d.``)
* rows without a gene symbol, empty rows and a header with a trailing space

Counts follow a negative binomial model, the standard model for RNA-seq
overdispersion, for 3 control vs 3 treated samples. The treated samples mimic
p53 activation: p53 target genes go up and cell-cycle genes go down.

Example:
    python generate_mock_rnaseq.py --output data/raw_counts.xlsx --seed 42
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

LOGGER = logging.getLogger("generate_mock_rnaseq")

GENE_COLUMN = "Gene_Symbol"
SAMPLES = ["Control_1", "Control_2", "Control_3", "Treated_1", "Treated_2", "Treated_3"]

# Real human (HGNC) gene symbols, grouped by how they respond to the treatment.
HOUSEKEEPING = """
ACTB GAPDH B2M RPLP0 TBP HPRT1 PPIA GUSB HMBS YWHAZ UBC TFRC PGK1 SDHA ALAS1
EEF1A1 RPL13A POLR2A
""".split()
P53_TARGETS = """
CDKN1A MDM2 BAX BBC3 PMAIP1 GADD45A FAS TP53I3 SESN1 SESN2 RRM2B ZMAT3
TNFRSF10B AEN PLK3 TIGAR DDB2 XPC FDXR RPS27L BTG2 CCNG1 TRIAP1 APAF1
""".split()
CELL_CYCLE = """
CCNB1 CCNB2 CDK1 PLK1 BUB1 AURKA AURKB CDC20 CCNA2 FOXM1 MKI67 TOP2A BIRC5
KIF20A CENPF E2F1 MCM2 CCNE1
""".split()
OTHER_GENES = """
TP53 MYC KRAS EGFR BRCA1 BRCA2 PTEN RB1 APC PIK3CA ERBB2 CDH1 VHL NRAS BRAF ATM
CHEK2 MDM4 NOTCH1 SMAD4 CTNNB1 IDH1 JAK2 ALK CD4 CD8A IL6 TNF STAT3 VEGFA HIF1A
NFKB1 IL1B CXCL8 TGFB1 MTOR AKT1 MAPK1 MAPK3 JUN FOS EGR1 HK2 LDHA PKM SLC2A1
G6PD IDH2
""".split()

TEXT_PLACEHOLDERS = ["NA", "n/a", "-", "n.d."]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate a messy, realistic RNA-seq raw-count Excel file\n"
        "for testing clean_data.py.",
        epilog="example:\n"
        "  python generate_mock_rnaseq.py --output data/raw_counts.xlsx",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("data/raw_counts.xlsx"),
        help="destination .xlsx file (default: %(default)s)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="random seed, for reproducible output (default: %(default)s)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="show debug messages"
    )
    args = parser.parse_args(argv)
    if args.output.suffix.lower() != ".xlsx":
        parser.error(f"--output must be an .xlsx file, got '{args.output}'")
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


def simulate_counts(rng: np.random.Generator) -> pd.DataFrame:
    """Simulate a clean gene x sample count matrix from a negative binomial model."""
    genes = HOUSEKEEPING + P53_TARGETS + CELL_CYCLE + OTHER_GENES
    n_variable = len(genes) - len(HOUSEKEEPING)

    # Mean expression: housekeeping genes are high and stable, the rest vary widely.
    base_mean = np.concatenate(
        [
            rng.lognormal(np.log(5000), 0.5, len(HOUSEKEEPING)),
            rng.lognormal(np.log(400), 1.2, n_variable),
        ]
    )
    # True treated-vs-control log2 fold changes.
    log2_fc = np.concatenate(
        [
            rng.normal(0.0, 0.1, len(HOUSEKEEPING)),
            rng.uniform(1.0, 3.0, len(P53_TARGETS)),
            -rng.uniform(1.0, 2.5, len(CELL_CYCLE)),
            rng.normal(0.0, 0.15, len(OTHER_GENES)),
        ]
    )
    is_treated = np.array([sample.startswith("Treated") for sample in SAMPLES])
    size_factors = rng.uniform(0.8, 1.25, len(SAMPLES))  # sequencing-depth differences

    mean = (
        base_mean[:, None]
        * size_factors
        * np.where(is_treated, 2.0 ** log2_fc[:, None], 1.0)
    )
    dispersion = 0.05 + 2.0 / mean  # DESeq2-like mean-dispersion trend
    size = 1.0 / dispersion
    counts = rng.negative_binomial(size, size / (size + mean))

    frame = pd.DataFrame(counts, columns=SAMPLES)
    frame.insert(0, GENE_COLUMN, genes)
    return frame


def make_messy(
    clean: pd.DataFrame, rng: np.random.Generator
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Inject the data-quality problems typically found in hand-edited spreadsheets."""
    issues: dict[str, int] = {}
    picks = rng.choice(len(clean), size=8, replace=False)

    # The same gene reported twice with different counts (e.g. two transcripts).
    second_feature = clean.iloc[picks[:5]].copy()
    second_feature[GENE_COLUMN] = second_feature[GENE_COLUMN].str.lower()
    second_feature[SAMPLES] = rng.poisson(0.3 * second_feature[SAMPLES].to_numpy())
    issues["Gene symbols repeated with different counts"] = len(second_feature)

    # Identical measurements pasted again with different formatting (" tp53").
    reformatted = clean.iloc[picks[5:]].copy()
    reformatted[GENE_COLUMN] = " " + reformatted[GENE_COLUMN].str.lower()
    issues["Re-formatted copies of existing rows"] = len(reformatted)

    messy = pd.concat([clean, second_feature, reformatted], ignore_index=True)
    messy[SAMPLES] = messy[SAMPLES].astype(object)  # cells may now hold text or NaN

    # Inconsistent capitalisation and stray whitespace in gene symbols.
    styles = [str.lower, str.capitalize, lambda s: f" {s}", lambda s: f"{s.lower()} "]
    restyled = rng.choice(len(clean), size=round(0.35 * len(clean)), replace=False)
    for row in restyled:
        style = styles[rng.integers(len(styles))]
        messy.loc[row, GENE_COLUMN] = style(messy.loc[row, GENE_COLUMN])
    issues["Mixed-case or space-padded gene symbols"] = len(restyled)

    # Blank cells and text placeholders in the count columns.
    n_cells = len(messy) * len(SAMPLES)
    n_blank = round(0.015 * n_cells)
    cells = rng.choice(n_cells, size=n_blank + len(TEXT_PLACEHOLDERS), replace=False)
    for i, cell in enumerate(cells):
        row, column = divmod(cell, len(SAMPLES))
        value = np.nan if i < n_blank else TEXT_PLACEHOLDERS[i - n_blank]
        messy.iat[row, 1 + column] = value
    issues["Blank count cells"] = n_blank
    issues["Text placeholders in count cells"] = len(TEXT_PLACEHOLDERS)

    # Rows without a gene symbol, plus completely empty rows.
    no_symbol = [[symbol, *rng.poisson(200, len(SAMPLES))] for symbol in (None, " ")]
    empty = [[None] * (len(SAMPLES) + 1)] * 2
    extra = pd.DataFrame(no_symbol + empty, columns=messy.columns, dtype=object)
    messy = pd.concat([messy, extra], ignore_index=True)
    issues["Rows without a gene symbol"] = len(no_symbol)
    issues["Completely empty rows"] = len(empty)

    # Exact duplicate rows, as if a block had been pasted twice.
    duplicates = messy.iloc[rng.choice(len(clean), size=8, replace=False)]
    messy = pd.concat([messy, duplicates], ignore_index=True)
    issues["Exact duplicate rows"] = len(duplicates)

    # Scatter the problems through the sheet and add a trailing space to a header.
    messy = messy.sample(frac=1.0, random_state=rng, ignore_index=True)
    messy = messy.rename(columns={"Treated_2": "Treated_2 "})
    issues["Column headers with trailing whitespace"] = 1
    return messy, issues


def main(argv: list[str] | None = None) -> int:
    """Run the generator and return a process exit code."""
    args = parse_args(argv)
    configure_logging(args.verbose)

    rng = np.random.default_rng(args.seed)
    clean = simulate_counts(rng)
    messy, issues = make_messy(clean, rng)

    try:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        messy.to_excel(args.output, index=False, sheet_name="raw_counts")
    except OSError as exc:
        LOGGER.error("Could not write %s: %s", args.output, exc)
        return 1

    LOGGER.info(
        "Simulated %d genes x %d samples (3 control vs 3 treated, seed=%d)",
        len(clean),
        len(SAMPLES),
        args.seed,
    )
    for issue, count in issues.items():
        LOGGER.info("  injected: %-44s %3d", issue, count)
    LOGGER.info("Wrote %d messy rows to %s", len(messy), args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
