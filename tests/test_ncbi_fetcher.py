"""Tests for Module B: PubMed record parsing and the CLI, fully offline.

NCBI is never contacted: the Bio.Entrez calls are replaced with fakes.
"""

from pathlib import Path
from urllib.error import URLError

import pandas as pd
import pytest
from Bio.Entrez.Parser import StringElement


def element(text: str, **attributes) -> StringElement:
    """A string with XML attributes, as returned by Bio.Entrez.read()."""
    return StringElement(text, tag="", attributes=attributes, key="")


@pytest.fixture
def record() -> dict:
    """A PubMed record shaped like Bio.Entrez.read() output."""
    return {
        "MedlineCitation": {
            "PMID": element("12345678"),
            "Article": {
                "ArticleTitle": element("The role of <i>TP53</i> in cancer."),
                "Abstract": {
                    "AbstractText": [
                        element("Why it matters.", Label="BACKGROUND"),
                        element("We  measured\n things.", Label="METHODS"),
                    ]
                },
                "AuthorList": [
                    {"LastName": "Smith", "Initials": "J"},
                    {"CollectiveName": "TP53 Consortium"},
                ],
                "Journal": {
                    "Title": "Nature genetics",
                    "JournalIssue": {"PubDate": {"MedlineDate": "2019 Jan-Feb"}},
                },
                "ELocationID": [element("10.1000/elocation", EIdType="doi")],
            },
        },
        "PubmedData": {
            "ArticleIdList": [
                element("12345678", IdType="pubmed"),
                element("10.1000/xyz", IdType="doi"),
            ]
        },
    }


def test_parse_article_flattens_record(fetch_abstracts, record):
    row = fetch_abstracts.parse_article(record)
    assert row == {
        "pmid": "12345678",
        "title": "The role of TP53 in cancer.",
        "abstract": "BACKGROUND: Why it matters. METHODS: We measured things.",
        "authors": "Smith J; TP53 Consortium",
        "journal": "Nature genetics",
        "year": "2019",
        "doi": "10.1000/xyz",
        "pubmed_url": "https://pubmed.ncbi.nlm.nih.gov/12345678/",
    }


def test_unlabelled_and_missing_abstracts(fetch_abstracts):
    unlabelled = {"Abstract": {"AbstractText": [element("Text.", Label="UNLABELLED")]}}
    assert fetch_abstracts.format_abstract(unlabelled) == "Text."
    assert fetch_abstracts.format_abstract({}) == ""


def test_doi_falls_back_to_electronic_location(fetch_abstracts, record):
    record["PubmedData"]["ArticleIdList"] = []
    assert fetch_abstracts.extract_doi(record) == "10.1000/elocation"


def test_default_output_is_named_after_the_query(fetch_abstracts):
    path = fetch_abstracts.default_output("TP53 AND Cancer")
    assert path == Path("output/pubmed_tp53_and_cancer.csv")


@pytest.mark.parametrize(
    "arguments",
    [
        ["-q", "TP53", "-e", "not-an-email"],
        ["-q", "TP53", "-e", "a@b.org", "--retmax", "0"],
        ["-q", "TP53", "-e", "a@b.org", "--retmax", "20000"],
    ],
)
def test_invalid_arguments_exit_with_usage_error(fetch_abstracts, arguments):
    with pytest.raises(SystemExit) as excinfo:
        fetch_abstracts.parse_args(arguments)
    assert excinfo.value.code == 2


class FakeHandle:
    """Stands in for the HTTP response returned by Entrez.esearch/efetch."""

    def __init__(self, kind: str):
        self.kind = kind

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        return False


def test_cli_writes_csv_without_network(fetch_abstracts, record, monkeypatch, tmp_path):
    entrez = fetch_abstracts.Entrez
    monkeypatch.setattr(entrez, "esearch", lambda **kwargs: FakeHandle("search"))
    monkeypatch.setattr(entrez, "efetch", lambda **kwargs: FakeHandle("fetch"))
    responses = {
        "search": {"Count": "1", "IdList": ["12345678"]},
        "fetch": {"PubmedArticle": [record]},
    }
    monkeypatch.setattr(entrez, "read", lambda handle: responses[handle.kind])
    output = tmp_path / "results.csv"

    exit_code = fetch_abstracts.main(["-q", "TP53", "-e", "a@b.org", "-o", str(output)])

    assert exit_code == 0
    table = pd.read_csv(output, dtype=str)
    assert table["pmid"].tolist() == ["12345678"]
    assert list(table.columns) == fetch_abstracts.COLUMNS


def test_cli_reports_network_errors(fetch_abstracts, monkeypatch, tmp_path):
    def offline(**kwargs):
        raise URLError("no internet")

    monkeypatch.setattr(fetch_abstracts.Entrez, "esearch", offline)
    exit_code = fetch_abstracts.main(
        ["-q", "TP53", "-e", "a@b.org", "-o", str(tmp_path / "x.csv")]
    )
    assert exit_code == 1
