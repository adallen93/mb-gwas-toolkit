"""Tests For The GWAS Toolkit."""

# conn = sqlite3.connect(":memory:")
# Parse_GWAS_Output("tests/data_for_testing.csv", ",", conn)
# gwas_obj = GWAS_Object(conn)
# print(gwas_obj.N_of_Significant_Tests)
# print(gwas_obj.Alpha_Level)
# print(gwas_obj.Manhattan_Plot)
# print(gwas_obj.QQ_Plot)

# Jonathan's Comments for Julia
# Show that the parse data function works well.
# Build tests to show all the GWAS_Object functions work.
#     The tricky thing: build tests proving that plot functions print plots.
#     This is pretty hard, but I think you can do it!
# Build tests that raises all error conditions.
# For some of these tests, you need to call data_for_testing.csv.
# But most of these tests don't need data_for_testing.csv.
# Instead, you can make your own toy data to prove, eg., that an error occurs.

import sqlite3
from typing import Any

import pytest
from gwas_toolkit import GWAS_Object, Parse_GWAS_Output


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


def test_qq_plot_displayed(caplog: Any) -> None:
    """Test to check if QQ plot is displayed."""
    conn = sqlite3.connect(":memory:")
    setup_database(conn)
    gwas = GWAS_Object(conn)

    gwas.Print_QQ_Plot()
    assert "A Separate Window Displaying he QQ Plot Was Opened." in caplog.text

    conn.close()


def test_manhattan_plot_displayed(caplog: Any) -> None:
    """Test to check if Manhattan plot is displayed."""
    conn = sqlite3.connect(":memory:")
    setup_database(conn)
    gwas = GWAS_Object(conn)

    gwas.Print_Manhattan_Plot()
    assert (
        "A Separate Window Displaying The Manhattan Plot Was Opened."
        in caplog.text
    )

    conn.close()
