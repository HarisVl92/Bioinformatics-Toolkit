"""Shared pytest fixtures.

The tool folders start with a digit (``1_data_wrangling`` ...), so they are not
importable Python packages. Each script is loaded straight from its file path.
"""

import importlib.util
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


def load_script(name: str, relative_path: str):
    """Import a stand-alone script as a module."""
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module  # dataclasses look the module up while importing
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def clean_data():
    return load_script("clean_data", "1_data_wrangling/clean_data.py")


@pytest.fixture(scope="session")
def generate_mock_rnaseq():
    return load_script(
        "generate_mock_rnaseq", "1_data_wrangling/generate_mock_rnaseq.py"
    )


@pytest.fixture(scope="session")
def fetch_abstracts():
    return load_script("fetch_abstracts", "2_ncbi_fetcher/fetch_abstracts.py")


@pytest.fixture(scope="session")
def plot_volcano():
    return load_script("plot_volcano", "3_data_visualization/plot_volcano.py")


@pytest.fixture(scope="session")
def generate_mock_de_results():
    return load_script(
        "generate_mock_de_results",
        "3_data_visualization/generate_mock_de_results.py",
    )
