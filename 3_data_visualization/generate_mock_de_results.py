#!/usr/bin/env python3
"""Generate mock differential-expression results in DESeq2's output format.

The table mimics ``results(dds)`` from DESeq2 for a treated-vs-control
comparison in which p53 has been activated: canonical p53 target genes are
up-regulated and cell-cycle genes (repressed through the p53-p21-DREAM
pathway) are down-regulated. All other genes are simulated background.
Like real DESeq2 output, some genes have no p-value (count outliers) or no
adjusted p-value (removed by independent filtering).

These are SIMULATED data for demonstrating plot_volcano.py, not the result
of a real experiment.

Example:
    python generate_mock_de_results.py --output data/mock_de_results.csv
"""

import argparse
import logging
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

LOGGER = logging.getLogger("generate_mock_de_results")

P53_TARGETS = """
CDKN1A MDM2 BAX BBC3 PMAIP1 GADD45A FAS TP53I3 SESN1 SESN2 RRM2B ZMAT3
TNFRSF10B AEN PLK3 TIGAR DDB2 XPC FDXR RPS27L BTG2 CCNG1 TRIAP1 APAF1
""".split()
CELL_CYCLE = """
CCNB1 CCNB2 CDK1 PLK1 BUB1 AURKA AURKB CDC20 CCNA2 FOXM1 MKI67 TOP2A BIRC5
KIF20A CENPF E2F1 MCM2 CCNE1
""".split()
MIN_GENES = 500


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate mock DESeq2-style differential-expression results\n"
        "for testing plot_volcano.py.",
        epilog="example:\n"
        "  python generate_mock_de_results.py --output data/mock_de_results.csv",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("data/mock_de_results.csv"),
        help="destination CSV file (default: %(default)s)",
    )
    parser.add_argument(
        "--n-genes",
        type=int,
        default=4000,
        help=f"number of genes to simulate, at least {MIN_GENES} "
        "(default: %(default)s)",
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
    if args.n_genes < MIN_GENES:
        parser.error(f"--n-genes must be at least {MIN_GENES}")
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


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR-adjusted p-values (NaN entries stay NaN)."""
    adjusted = np.full(pvalues.shape, np.nan)
    observed = ~np.isnan(pvalues)
    values = pvalues[observed]
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]  # enforce monotonicity
    result = np.empty(len(values))
    result[order] = np.clip(ranked, 0.0, 1.0)
    adjusted[observed] = result
    return adjusted


def simulate_results(n_genes: int, rng: np.random.Generator) -> pd.DataFrame:
    """Simulate a DESeq2-like results table for a p53-activation experiment."""
    n_known = len(P53_TARGETS) + len(CELL_CYCLE)
    n_background = n_genes - n_known
    genes = (
        P53_TARGETS + CELL_CYCLE + [f"GENE{i:05d}" for i in range(1, n_background + 1)]
    )

    # Background genes: ~90% unchanged, ~10% moderately up or down.
    background_lfc = np.zeros(n_background)
    regulated = rng.random(n_background) < 0.10
    background_lfc[regulated] = rng.choice([-1.0, 1.0], regulated.sum()) * rng.uniform(
        0.5, 2.5, regulated.sum()
    )
    true_lfc = np.concatenate(
        [
            rng.uniform(1.5, 4.0, len(P53_TARGETS)),
            -rng.uniform(1.2, 3.0, len(CELL_CYCLE)),
            background_lfc,
        ]
    )
    base_mean = np.concatenate(
        [
            rng.lognormal(np.log(2000), 0.6, n_known),  # well-expressed known genes
            rng.lognormal(np.log(300), 1.6, n_background),
        ]
    )

    # Wald test: low-count genes have noisier fold-change estimates.
    lfc_se = 0.25 + 2.5 / np.sqrt(base_mean)
    log2_fc = true_lfc + rng.normal(0.0, lfc_se)
    stat = log2_fc / lfc_se
    pvalue = np.vectorize(math.erfc)(np.abs(stat) / math.sqrt(2.0))  # two-sided

    # DESeq2 sets p-values of count outliers to NA ...
    outlier = rng.random(n_genes) < 0.003
    stat[outlier] = np.nan
    pvalue[outlier] = np.nan
    # ... and excludes low-count genes from the FDR correction (padj = NA).
    filtered = base_mean < np.quantile(base_mean, 0.08)
    padj = benjamini_hochberg(np.where(filtered, np.nan, pvalue))

    results = pd.DataFrame(
        {
            "gene": genes,
            "baseMean": base_mean,
            "Log2FoldChange": log2_fc,
            "lfcSE": lfc_se,
            "stat": stat,
            "pvalue": pvalue,
            "padj": padj,
        }
    )
    return results.sort_values("pvalue", ignore_index=True)  # NA p-values last


def main(argv: list[str] | None = None) -> int:
    """Run the generator and return a process exit code."""
    args = parse_args(argv)
    configure_logging(args.verbose)

    results = simulate_results(args.n_genes, np.random.default_rng(args.seed))
    try:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        results.to_csv(args.output, index=False, float_format="%.6g")
    except OSError as exc:
        LOGGER.error("Could not write %s: %s", args.output, exc)
        return 1

    significant = int((results["padj"] < 0.05).sum())
    LOGGER.info(
        "Simulated %d genes (seed=%d): %d p53 targets up, %d cell-cycle genes down",
        len(results),
        args.seed,
        len(P53_TARGETS),
        len(CELL_CYCLE),
    )
    LOGGER.info(
        "%d genes with padj < 0.05; %d without p-value, %d without padj (as in DESeq2)",
        significant,
        results["pvalue"].isna().sum(),
        results["padj"].isna().sum(),
    )
    LOGGER.info("Wrote mock DESeq2 results to %s", args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
