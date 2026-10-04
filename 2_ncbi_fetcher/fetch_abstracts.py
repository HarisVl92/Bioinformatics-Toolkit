#!/usr/bin/env python3
"""Search PubMed and export article metadata and abstracts to CSV.

Uses the NCBI E-utilities through Biopython's ``Bio.Entrez``: ESearch finds the
PMIDs that match a query, EFetch downloads the full records in batches, and
each record is flattened into one CSV row (PMID, title, abstract, authors,
journal, year, DOI and a link to the PubMed page).

NCBI asks every user to identify themselves with an e-mail address and limits
traffic to 3 requests per second (10 with an API key). Biopython enforces the
rate limit and retries temporary server errors automatically.

Example:
    python fetch_abstracts.py --query "TP53 AND Cancer" --retmax 15 --email you@uni.edu
"""

import argparse
import logging
import os
import re
import sys
from pathlib import Path
from urllib.error import HTTPError, URLError

import pandas as pd
from Bio import Entrez

LOGGER = logging.getLogger("fetch_abstracts")

TOOL_NAME = "bioinformatics-toolkit"
BATCH_SIZE = 100  # records per EFetch request
MAX_RECORDS = 10_000  # ESearch cannot page beyond 10,000 PubMed records
COLUMNS = [
    "pmid",
    "title",
    "abstract",
    "authors",
    "journal",
    "year",
    "doi",
    "pubmed_url",
]
EMAIL_PATTERN = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")
YEAR_PATTERN = re.compile(r"\b(?:19|20)\d{2}\b")
# Inline formatting that PubMed keeps in titles/abstracts, e.g. <i>TP53</i>, <sup>2</sup>
INLINE_TAG_PATTERN = re.compile(r"</?(?:i|b|u|em|strong|sup|sub)(?:\s[^>]*)?>", re.I)


def positive_int(value: str) -> int:
    """argparse type: an integer >= 1."""
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError(f"must be a positive integer, got {value}")
    return number


def email_address(value: str) -> str:
    """argparse type: a syntactically valid e-mail address."""
    if not EMAIL_PATTERN.match(value):
        raise argparse.ArgumentTypeError(f"'{value}' is not a valid e-mail address")
    return value


def default_output(query: str) -> Path:
    """Build a file name from the query, e.g. 'TP53 AND Cancer' -> pubmed_tp53_and_cancer.csv."""
    slug = re.sub(r"[^a-z0-9]+", "_", query.lower()).strip("_")[:60] or "results"
    return Path("output") / f"pubmed_{slug}.csv"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse and validate command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Search PubMed and save article metadata and abstracts to CSV.",
        epilog="examples:\n"
        '  python fetch_abstracts.py --query "TP53 AND Cancer" --retmax 15 '
        "--email you@uni.edu\n"
        '  python fetch_abstracts.py -q "BRCA1[tiab] AND 2020:2025[dp]" -n 200 '
        "-e you@uni.edu --sort pub_date",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-q",
        "--query",
        required=True,
        help='PubMed search term, e.g. "TP53 AND Cancer" (full PubMed syntax, '
        'such as "BRCA1[tiab] AND 2020:2025[dp]")',
    )
    parser.add_argument(
        "-n",
        "--retmax",
        type=positive_int,
        default=20,
        help="maximum number of records to download, up to 10,000 "
        "(default: %(default)s)",
    )
    parser.add_argument(
        "-e",
        "--email",
        type=email_address,
        required=True,
        help="your e-mail address (required by NCBI's E-utilities usage policy)",
    )
    parser.add_argument(
        "--api-key",
        help="optional NCBI API key; raises the rate limit from 3 to 10 requests/s "
        "(default: $NCBI_API_KEY if set)",
    )
    parser.add_argument(
        "--sort",
        choices=("relevance", "pub_date"),
        default="relevance",
        help="'relevance' (PubMed Best Match) or 'pub_date' (newest first) "
        "(default: %(default)s)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="destination CSV (default: output/pubmed_<query>.csv)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="show debug messages"
    )
    args = parser.parse_args(argv)
    if args.retmax > MAX_RECORDS:
        parser.error(f"--retmax cannot exceed {MAX_RECORDS:,} (an E-utilities limit)")
    if args.output is None:
        args.output = default_output(args.query)
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


def search_pubmed(query: str, retmax: int, sort: str) -> tuple[int, list[str]]:
    """Run ESearch; return the total number of hits and the first PMIDs."""
    with Entrez.esearch(db="pubmed", term=query, retmax=retmax, sort=sort) as handle:
        result = Entrez.read(handle)
    LOGGER.debug("PubMed query translation: %s", result.get("QueryTranslation", ""))
    return int(result["Count"]), list(result["IdList"])


def fetch_records(pmids: list[str]) -> list:
    """Download full PubMed records with EFetch, in batches."""
    records = []
    for start in range(0, len(pmids), BATCH_SIZE):
        batch = pmids[start : start + BATCH_SIZE]
        LOGGER.info(
            "Downloading records %d-%d of %d",
            start + 1,
            start + len(batch),
            len(pmids),
        )
        with Entrez.efetch(db="pubmed", id=",".join(batch), retmode="xml") as handle:
            records.extend(Entrez.read(handle).get("PubmedArticle", []))
    return records


def clean_text(value: object) -> str:
    """Remove inline HTML formatting and collapse whitespace into single spaces."""
    text = INLINE_TAG_PATTERN.sub("", str(value))
    return re.sub(r"\s+", " ", text).strip()


def format_abstract(article: dict) -> str:
    """Join abstract sections, keeping the labels of structured abstracts."""
    parts = []
    for section in article.get("Abstract", {}).get("AbstractText", []):
        label = getattr(section, "attributes", {}).get("Label", "")
        text = clean_text(section)
        if not text:
            continue
        # NLM marks sections without a heading as "UNLABELLED"; omit that label.
        has_label = label and label.upper() != "UNLABELLED"
        parts.append(f"{label}: {text}" if has_label else text)
    return " ".join(parts)


def format_authors(authors: list) -> str:
    """Format authors as 'Lastname Initials; ...' (consortia by their name)."""
    names = []
    for author in authors:
        if "LastName" in author:
            names.append(f"{author['LastName']} {author.get('Initials', '')}".strip())
        elif "CollectiveName" in author:
            names.append(clean_text(author["CollectiveName"]))
    return "; ".join(names)


def extract_year(article: dict) -> str:
    """Publication year from the journal issue date, or the electronic date."""
    pub_date = article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
    candidates = [pub_date.get("Year", ""), pub_date.get("MedlineDate", "")]
    candidates += [date.get("Year", "") for date in article.get("ArticleDate", [])]
    for text in candidates:
        match = YEAR_PATTERN.search(str(text))  # MedlineDate looks like "2019 Jan-Feb"
        if match:
            return match.group(0)
    return ""


def extract_doi(record: dict) -> str:
    """DOI from the PubMed article IDs, falling back to the electronic location."""
    for article_id in record.get("PubmedData", {}).get("ArticleIdList", []):
        if getattr(article_id, "attributes", {}).get("IdType") == "doi":
            return str(article_id)
    article = record["MedlineCitation"]["Article"]
    for location in article.get("ELocationID", []):
        if getattr(location, "attributes", {}).get("EIdType") == "doi":
            return str(location)
    return ""


def parse_article(record: dict) -> dict[str, str]:
    """Flatten one PubmedArticle record into a CSV row."""
    citation = record["MedlineCitation"]
    article = citation["Article"]
    pmid = str(citation["PMID"])
    return {
        "pmid": pmid,
        "title": clean_text(article.get("ArticleTitle", "")),
        "abstract": format_abstract(article),
        "authors": format_authors(article.get("AuthorList", [])),
        "journal": clean_text(article.get("Journal", {}).get("Title", "")),
        "year": extract_year(article),
        "doi": extract_doi(record),
        "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
    }


def log_preview(table: pd.DataFrame, n_rows: int = 3) -> None:
    """Show the first few titles so the user can sanity-check the results."""
    for row in table.head(n_rows).itertuples():
        title = row.title if len(row.title) <= 80 else row.title[:77] + "..."
        LOGGER.info("  [%s] %s (PMID %s)", row.year or "n.d.", title, row.pmid)


def main(argv: list[str] | None = None) -> int:
    """Search PubMed, download the records and save them to CSV."""
    args = parse_args(argv)
    configure_logging(args.verbose)
    Entrez.email = args.email
    Entrez.tool = TOOL_NAME
    Entrez.api_key = args.api_key or os.environ.get("NCBI_API_KEY")

    try:
        total, pmids = search_pubmed(args.query, args.retmax, args.sort)
        if not pmids:
            LOGGER.warning("No PubMed records match %r", args.query)
        else:
            LOGGER.info(
                "PubMed: %s records match %r; downloading %d",
                f"{total:,}",
                args.query,
                len(pmids),
            )
        records = fetch_records(pmids)
    except HTTPError as exc:
        LOGGER.error(
            "NCBI returned HTTP %s (%s). Try again later.", exc.code, exc.reason
        )
        return 1
    except (URLError, TimeoutError, ConnectionError) as exc:
        reason = getattr(exc, "reason", exc)
        LOGGER.error("Could not reach NCBI (%s). Check your connection.", reason)
        return 1
    except RuntimeError as exc:  # error messages that NCBI embeds in its XML
        LOGGER.error("NCBI reported an error: %s", exc)
        return 1
    except ValueError as exc:  # malformed or unexpected XML
        LOGGER.error("Unexpected response from NCBI: %s", exc)
        return 1
    except KeyboardInterrupt:
        LOGGER.error("Interrupted by user")
        return 130

    table = pd.DataFrame([parse_article(record) for record in records], columns=COLUMNS)
    try:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(args.output, index=False)
    except OSError as exc:
        LOGGER.error("Could not write %s: %s", args.output, exc)
        return 1

    with_abstract = int((table["abstract"] != "").sum())
    LOGGER.info(
        "Saved %d records (%d with an abstract) to %s",
        len(table),
        with_abstract,
        args.output,
    )
    log_preview(table)
    return 0


if __name__ == "__main__":
    sys.exit(main())
