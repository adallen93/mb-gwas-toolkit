"""Tests For The GWAS Toolkit."""

import sqlite3
import warnings
from io import BytesIO
from typing import Any

import matplotlib.pyplot as plt
import pytest
from gwas_toolkit import GWAS_Object, Parse_GWAS_Output


@pytest.fixture(autouse=True)
def no_show_warnings() -> Any:
    """Ignore Warnings."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        yield


def test_parse_data() -> None:
    """A Test To Determine If Parse Data Interacts With SQL As Intended."""
    conn = sqlite3.connect(":memory:")
    Parse_GWAS_Output("tests/test_data.csv", ",", conn)
    c = conn.cursor()
    c.execute("SELECT * FROM gwas")
    allelle = c.fetchall()
    assert len(allelle) > 0
    conn.close()


def setup_database(conn: sqlite3.Connection) -> None:
    """Populates the SQLite database with data from the test_data.csv file."""
    c = conn.cursor()
    c.execute("""CREATE TABLE IF NOT EXISTS gwas (
                 MarkerID TEXT,
                 Chromosome TEXT,
                 Location REAL,
                 PValue REAL)""")
    c.executemany(
        """INSERT INTO gwas
        (MarkerID, Chromosome, Location, PValue)
        VALUES (?, ?, ?, ?)""",
        [
            ("Marker1", "1", "100", "0.1"),
            ("Marker2", "2", "200", "0.001"),
        ],
    )
    conn.commit()


def test_alpha_level() -> None:
    """Test to check the alpha level calculation."""
    conn = sqlite3.connect(":memory:")
    setup_database(conn)
    gwas = GWAS_Object(conn)
    assert gwas.Alpha_Level == 0.05
    conn.close()


def test_qq_plot_displayed() -> None:
    """Test to check if QQ plot is displayed."""
    conn = sqlite3.connect(":memory:")
    setup_database(conn)
    gwas = GWAS_Object(conn)
    buffer = BytesIO()
    plt.switch_backend("agg")
    plt.figure()
    gwas.Print_QQ_Plot()
    plt.savefig(buffer, format="png")
    plt.close()
    assert buffer.getbuffer().nbytes > 0


def test_manhattan_plot_displayed() -> None:
    """Test to check if Manhattan plot is displayed."""
    conn = sqlite3.connect(":memory:")
    setup_database(conn)
    gwas = GWAS_Object(conn)
    buffer = BytesIO()
    plt.switch_backend("agg")
    plt.figure()
    gwas.Print_Manhattan_Plot()
    plt.savefig(buffer, format="png")
    plt.close()
    assert buffer.getbuffer().nbytes > 0
