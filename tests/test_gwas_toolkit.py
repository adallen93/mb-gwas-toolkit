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


@pytest.fixture
def gwas_object() -> Any:
    """Fixture to create a GWAS_Object instance for testing."""
    conn: sqlite3.Connection = sqlite3.connect(":memory:")
    Parse_GWAS_Output("tests/test_data.csv", ",", conn)
    gwas: GWAS_Object = GWAS_Object(conn)
    yield gwas
    conn.close()


def test_alpha_level(gwas_object: GWAS_Object) -> None:
    """Test to check the alpha level calculation."""
    # Assert that the alpha level is 0.05
    assert gwas_object.Alpha_Level == 0.05


def test_qq_plot_displayed(gwas_object: GWAS_Object, caplog: Any) -> None:
    """Test to check if QQ plot is displayed."""
    # Assert that calling Print_QQ_Plot displays the QQ plot
    gwas_object.Print_QQ_Plot()
    assert "A Separate Window Displaying he QQ Plot Was Opened." in caplog.text


def test_manhattan_plot_displayed(
    gwas_object: GWAS_Object, caplog: Any
) -> None:
    """Test to check if Manhattan plot is displayed."""
    # Assert that calling Print_Manhattan_Plot displays the Manhattan plot
    gwas_object.Print_Manhattan_Plot()
    assert (
        "A Separate Window Displaying The Manhattan Plot Was Opened."
        in caplog.text
    )
