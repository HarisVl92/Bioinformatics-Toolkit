#!/usr/bin/env bash
# Reproduce every demo output in this repository with a single command.
#
# Usage (with the virtual environment activated):
#   NCBI_EMAIL=you@university.edu bash run_demo.sh
#
# Set PYTHON=/path/to/python to use a specific interpreter instead.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python}"
: "${NCBI_EMAIL:?set NCBI_EMAIL to your e-mail address, e.g. NCBI_EMAIL=you@uni.edu bash run_demo.sh}"

BOLD="" RESET=""
if [[ -t 1 ]]; then  # bold step headers only when writing to a terminal
    BOLD="$(tput bold 2>/dev/null || true)"
    RESET="$(tput sgr0 2>/dev/null || true)"
fi
step() { printf '\n%s==> %s%s\n' "$BOLD" "$1" "$RESET"; }

step "[1/3] Data wrangling: messy Excel sheet -> clean count matrix"
cd "$ROOT/1_data_wrangling"
"$PYTHON" generate_mock_rnaseq.py --output data/raw_counts.xlsx
"$PYTHON" clean_data.py --input data/raw_counts.xlsx --output output/clean_counts.csv

step "[2/3] PubMed fetcher: \"TP53 AND Cancer\" -> CSV of abstracts"
cd "$ROOT/2_ncbi_fetcher"
"$PYTHON" fetch_abstracts.py --query "TP53 AND Cancer" --retmax 15 --email "$NCBI_EMAIL"

step "[3/3] Visualization: differential-expression results -> volcano plot"
cd "$ROOT/3_data_visualization"
"$PYTHON" generate_mock_de_results.py --output data/mock_de_results.csv
"$PYTHON" plot_volcano.py --input data/mock_de_results.csv \
    --output output/volcano_plot.png --title "Treated vs Control (simulated data)"

step "Done: results are in 1_data_wrangling/output, 2_ncbi_fetcher/output and 3_data_visualization/output"
