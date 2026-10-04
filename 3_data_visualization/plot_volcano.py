#!/usr/bin/env python3
"""Create a publication-ready volcano plot from differential-expression results.

Reads a table with a log2 fold-change column and a p-value column (DESeq2,
edgeR, limma or similar), classifies every gene as up-regulated,
down-regulated or not significant, labels the top hits and saves a
high-resolution figure (PNG, PDF, SVG or TIFF).

Column names are matched case-insensitively, so DESeq2's "log2FoldChange"
works with the default --fc-col.

Example:
    python plot_volcano.py --input data/mock_de_results.csv --output output/volcano_plot.png
"""

import argparse
import logging
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")  # render without a display (servers, CI, SSH sessions)

import matplotlib.pyplot as plt  # noqa: E402  (must follow matplotlib.use)
import seaborn as sns  # noqa: E402
from adjustText import adjust_text  # noqa: E402
from matplotlib import patheffects  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

LOGGER = logging.getLogger("plot_volcano")

UP, DOWN, NOT_SIGNIFICANT = "Up-regulated", "Down-regulated", "Not significant"
GROUP_ORDER = [UP, DOWN, NOT_SIGNIFICANT]
# Colour-blind-safe pair, validated under protanopia/deuteranopia simulation
# (OKLab delta-E >= 21) with >= 3:1 contrast on white; grey lets the bulk recede.
PALETTE = {UP: "#e34948", DOWN: "#2a78d6", NOT_SIGNIFICANT: "#c3c2b7"}
MARKER_SIZES = {UP: 22, DOWN: 22, NOT_SIGNIFICANT: 10}
INK, SECONDARY_INK, MUTED = "#0b0b0b", "#52514e", "#898781"
GENE_COLUMN_CANDIDATES = ("gene", "gene_symbol", "symbol", "gene_name", "gene_id")
IMAGE_SUFFIXES = {".png", ".pdf", ".svg", ".tif", ".tiff"}


def positive_float(value: str) -> float:
    """argparse type: a number > 0."""
    number = float(value)
    if not number > 0:  # also rejects NaN
        raise argparse.ArgumentTypeError(f"must be a positive number, got {value}")
    return number


def probability(value: str) -> float:
    """argparse type: a number strictly between 0 and 1."""
    number = float(value)
    if not 0 < number < 1:
        raise argparse.ArgumentTypeError(f"must be between 0 and 1, got {value}")
    return number


def non_negative_int(value: str) -> int:
    """argparse type: an integer >= 0."""
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError(f"must be 0 or more, got {value}")
    return number


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse and validate command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Create a publication-ready volcano plot from differential-\n"
        "expression results (DESeq2, edgeR, limma, ...).",
        epilog="examples:\n"
        "  python plot_volcano.py --input data/mock_de_results.csv "
        "--output output/volcano_plot.png\n"
        "  python plot_volcano.py -i deseq2_results.csv -o volcano.pdf "
        "--pval-col padj --top-n 20",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="CSV (or .tsv) with one row per gene",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("volcano_plot.png"),
        help="output image; the format follows the extension: .png, .pdf, .svg, "
        ".tif (default: %(default)s)",
    )
    parser.add_argument(
        "--fc-col",
        default="Log2FoldChange",
        help="log2 fold-change column, case-insensitive (default: %(default)s)",
    )
    parser.add_argument(
        "--pval-col",
        default="pvalue",
        help="p-value column; use 'padj' to call significance on FDR-adjusted "
        "p-values (default: %(default)s)",
    )
    parser.add_argument(
        "--gene-col",
        help="column with gene names for the labels (default: auto-detect)",
    )
    parser.add_argument(
        "--fc-threshold",
        type=positive_float,
        default=1.0,
        help="minimum |log2 fold change|; 1.0 means a two-fold change "
        "(default: %(default)s)",
    )
    parser.add_argument(
        "--pval-threshold",
        type=probability,
        default=0.05,
        help="significance cut-off (default: %(default)s)",
    )
    parser.add_argument(
        "--top-n",
        type=non_negative_int,
        default=10,
        help="number of most significant genes to label, split between up and "
        "down; 0 disables labels (default: %(default)s)",
    )
    parser.add_argument(
        "--title", default="Volcano plot", help="figure title (default: %(default)s)"
    )
    parser.add_argument(
        "--dpi",
        type=non_negative_int,
        default=300,
        help="resolution for PNG/TIFF output (default: %(default)s)",
    )
    parser.add_argument(
        "--figsize",
        type=positive_float,
        nargs=2,
        default=(7.0, 6.0),
        metavar=("WIDTH", "HEIGHT"),
        help="figure size in inches (default: 7 6)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="show debug messages"
    )
    args = parser.parse_args(argv)
    if args.output.suffix.lower() not in IMAGE_SUFFIXES:
        parser.error(
            f"unsupported image format '{args.output.suffix}'; "
            f"use one of {', '.join(sorted(IMAGE_SUFFIXES))}"
        )
    if args.dpi < 50:
        parser.error("--dpi must be at least 50")
    return args


def configure_logging(verbose: bool) -> None:
    """Log to stderr as 'time | level | message'.

    Only this tool logs at INFO/DEBUG; third-party libraries report warnings only
    (otherwise e.g. fontTools floods the console when saving PDFs).
    """
    logging.basicConfig(
        level=logging.WARNING,
        format="%(asctime)s | %(levelname)-8s | %(message)s",
        datefmt="%H:%M:%S",
    )
    LOGGER.setLevel(logging.DEBUG if verbose else logging.INFO)


def load_table(path: Path) -> pd.DataFrame:
    """Read a CSV (or a tab-separated .tsv/.txt) file."""
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path}")
    separator = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=separator)


def match_column(columns: pd.Index, name: str) -> str | None:
    """Return the column called ``name`` ignoring case and spaces, or None."""
    wanted = name.strip().lower()
    return next((c for c in columns if str(c).strip().lower() == wanted), None)


def find_column(frame: pd.DataFrame, name: str) -> str:
    """Like match_column, but raise a helpful error when the column is missing."""
    column = match_column(frame.columns, name)
    if column is None:
        available = ", ".join(map(str, frame.columns))
        raise ValueError(f"column '{name}' not found; available columns: {available}")
    return column


def find_gene_column(frame: pd.DataFrame, requested: str | None) -> str | None:
    """Pick the gene-label column: the requested one, a common name, or the first
    text column (e.g. the unnamed row-name column written by R's write.csv)."""
    if requested:
        return find_column(frame, requested)
    for candidate in GENE_COLUMN_CANDIDATES:
        column = match_column(frame.columns, candidate)
        if column is not None:
            return column
    first = frame.columns[0]
    return None if pd.api.types.is_numeric_dtype(frame[first]) else first


def prepare_data(
    frame: pd.DataFrame,
    fc_col: str,
    pval_col: str,
    gene_col: str | None,
    fc_threshold: float,
    pval_threshold: float,
) -> pd.DataFrame:
    """Validate the inputs and classify each gene as up, down or not significant."""
    data = pd.DataFrame(
        {
            "gene": frame[gene_col].astype(str) if gene_col else "",
            "log2fc": pd.to_numeric(frame[fc_col], errors="coerce"),
            "pvalue": pd.to_numeric(frame[pval_col], errors="coerce"),
        }
    )
    missing = data[["log2fc", "pvalue"]].isna().any(axis="columns")
    if missing.any():
        LOGGER.warning(
            "Skipping %d gene(s) without a fold change or p-value "
            "(e.g. DESeq2 outliers or filtered genes)",
            missing.sum(),
        )
        data = data.loc[~missing]
    if data.empty:
        raise ValueError("no genes have both a fold change and a p-value")
    if not data["pvalue"].between(0, 1).all():
        raise ValueError(
            f"column '{pval_col}' has values outside 0-1; is it a p-value column?"
        )

    zero = data["pvalue"] == 0
    if zero.any():  # DESeq2 reports p = 0 when the true value underflows
        floor = data.loc[~zero, "pvalue"].min() if (~zero).any() else 1e-300
        LOGGER.warning(
            "%d p-value(s) are exactly 0; plotting them at the smallest "
            "non-zero p-value (%.1e)",
            zero.sum(),
            floor,
        )
        data = data.assign(pvalue=data["pvalue"].mask(zero, floor))

    significant = data["pvalue"] < pval_threshold
    regulation = np.select(
        [
            significant & (data["log2fc"] >= fc_threshold),
            significant & (data["log2fc"] <= -fc_threshold),
        ],
        [UP, DOWN],
        default=NOT_SIGNIFICANT,
    )
    return data.assign(neg_log10_p=-np.log10(data["pvalue"]), regulation=regulation)


def pick_labels(data: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """The most significant genes to label, split evenly between up and down."""
    up = data[data["regulation"] == UP].nsmallest(math.ceil(top_n / 2), "pvalue")
    down = data[data["regulation"] == DOWN].nsmallest(top_n // 2, "pvalue")
    return pd.concat([up, down])


def draw_volcano(
    data: pd.DataFrame,
    *,
    fc_threshold: float,
    pval_threshold: float,
    p_label: str,
    title: str,
    top_n: int,
    figsize: tuple[float, float],
) -> plt.Figure:
    """Draw the volcano plot and return the figure."""
    sns.set_theme(context="paper", style="ticks", font_scale=1.2)
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,  # embed TrueType fonts: editable text, journal-safe
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig, ax = plt.subplots(figsize=figsize)

    counts = data["regulation"].value_counts()
    names = {group: f"{group} (n={counts.get(group, 0):,})" for group in GROUP_ORDER}
    # Grey points first, so the coloured (significant) points are drawn on top.
    plot_data = data.assign(group=data["regulation"].map(names)).sort_values(
        "regulation", key=lambda column: column.ne(NOT_SIGNIFICANT), kind="stable"
    )
    order = [names[group] for group in GROUP_ORDER]
    sns.scatterplot(
        data=plot_data,
        x="log2fc",
        y="neg_log10_p",
        hue="group",
        hue_order=order,
        palette={names[group]: PALETTE[group] for group in GROUP_ORDER},
        size="group",
        size_order=order,
        sizes={names[group]: MARKER_SIZES[group] for group in GROUP_ORDER},
        alpha=0.85,
        linewidth=0.3,
        edgecolor="white",
        rasterized=True,  # keeps PDF/SVG small while text stays vector
        legend=False,
        ax=ax,
    )

    # Threshold guides behind the points.
    for x in (-fc_threshold, fc_threshold):
        ax.axvline(x, color=MUTED, linestyle="--", linewidth=0.8, zorder=0)
    ax.axhline(
        -math.log10(pval_threshold),
        color=MUTED,
        linestyle="--",
        linewidth=0.8,
        zorder=0,
    )

    # Symmetric x-axis so up- and down-regulation are visually comparable.
    limit = max(data["log2fc"].abs().max(), fc_threshold) * 1.1
    ax.set_xlim(-limit, limit)
    ax.set_ylim(0, data["neg_log10_p"].max() * 1.1)
    ax.set_xlabel("log₂ fold change")
    ax.set_ylabel(f"−log₁₀({p_label})")
    ax.set_title(title, loc="left", fontweight="bold", color=INK, pad=22)
    ax.text(
        0,
        1.02,
        f"Significant: |log₂ FC| ≥ {fc_threshold:g} and {p_label} < "
        f"{pval_threshold:g}",
        transform=ax.transAxes,
        fontsize=9,
        color=SECONDARY_INK,
    )
    # Hand-built legend: equal-size markers are easier to read than the
    # point sizes used in the plot.
    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            markersize=6,
            markerfacecolor=PALETTE[group],
            markeredgewidth=0,
            label=names[group],
        )
        for group in GROUP_ORDER
    ]
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=3,
        frameon=False,
        fontsize=9,
        handletextpad=0.1,
        columnspacing=1.0,
    )
    sns.despine(ax=ax)

    labels = pick_labels(data, top_n)
    if not labels.empty:
        # A white halo keeps labels legible where they cross points or leader lines.
        halo = [patheffects.withStroke(linewidth=3, foreground="white")]
        texts = [
            ax.text(
                row.log2fc,
                row.neg_log10_p,
                row.gene,
                fontsize=8,
                color=INK,
                zorder=5,
                path_effects=halo,
            )
            for row in labels.itertuples()
        ]
        adjust_text(
            texts,
            x=data["log2fc"].to_numpy(),
            y=data["neg_log10_p"].to_numpy(),
            ax=ax,
            expand=(1.3, 1.6),  # extra padding around labels
            force_text=(0.3, 0.6),  # push overlapping labels further apart
            arrowprops={"arrowstyle": "-", "color": MUTED, "lw": 0.6},
        )
    return fig


def save_figure(fig: plt.Figure, path: Path, dpi: int) -> None:
    """Save the figure (format from the file extension), creating folders."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    """Build the volcano plot and return a process exit code."""
    args = parse_args(argv)
    configure_logging(args.verbose)

    try:
        table = load_table(args.input)
        fc_col = find_column(table, args.fc_col)
        pval_col = find_column(table, args.pval_col)
        gene_col = find_gene_column(table, args.gene_col)
        LOGGER.info(
            "Loaded %d genes from %s (fold change: '%s', p-value: '%s', labels: '%s')",
            len(table),
            args.input,
            fc_col,
            pval_col,
            gene_col,
        )
        if gene_col is None and args.top_n:
            LOGGER.warning("No gene-name column found; use --gene-col to add labels")
        data = prepare_data(
            table, fc_col, pval_col, gene_col, args.fc_threshold, args.pval_threshold
        )
        is_adjusted = any(key in pval_col.lower() for key in ("adj", "fdr", "qval"))
        fig = draw_volcano(
            data,
            fc_threshold=args.fc_threshold,
            pval_threshold=args.pval_threshold,
            p_label="adjusted p-value" if is_adjusted else "p-value",
            title=args.title,
            top_n=args.top_n if gene_col else 0,
            figsize=tuple(args.figsize),
        )
        save_figure(fig, args.output, args.dpi)
    except FileNotFoundError as exc:
        LOGGER.error("%s", exc)
        return 1
    except ValueError as exc:
        LOGGER.error("Invalid input: %s", exc)
        return 1
    except OSError as exc:
        LOGGER.error("File error: %s", exc)
        return 1
    except KeyboardInterrupt:
        LOGGER.error("Interrupted by user")
        return 130

    counts = data["regulation"].value_counts()
    LOGGER.info(
        "Up-regulated: %d | Down-regulated: %d | Not significant: %d",
        counts.get(UP, 0),
        counts.get(DOWN, 0),
        counts.get(NOT_SIGNIFICANT, 0),
    )
    LOGGER.info("Saved volcano plot (%d dpi) to %s", args.dpi, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
