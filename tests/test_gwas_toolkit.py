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

import pytest
from gwas_toolkit import (
    Benjamini_Hochberg_Procedure,
    DataNotFoundError,
    GWAS_Object,
    Parse_GWAS_Output,
)


def test_print_manhattan_plot() -> None:
    """Tests the Print_QQ_Plot method of the GWAS class."""
    # First create a sqlite3 connection with memory
    conn = sqlite3.connect(":memory:")
    c = conn.cursor()

    # Replace any existing table with an emtpy table
    c.execute("""CREATE TABLE IF NOT EXISTS gwas (
              MarkerID TEXT,
              Chromosome TEXT,
              Location REAL,
              PValue REAL
              )""")

    # Create a GWAS_Object tied to the empty data
    gwas = GWAS_Object(conn)

    # Test that this method raises an error
    with pytest.raises(DataNotFoundError):
        gwas.Print_Manhattan_Plot()

    # Because of the nature of this method, we are unsure of how to
    # test the remaining parts

    # Close the connection
    conn.close()


def test_print_qq_plot() -> None:
    """Tests the Print_QQ_Plot method of the GWAS class."""
    # First create a sqlite3 connection with memory
    conn = sqlite3.connect(":memory:")
    c = conn.cursor()

    # Replace any existing table with an emtpy table
    c.execute("""CREATE TABLE IF NOT EXISTS gwas (
              MarkerID TEXT,
              Chromosome TEXT,
              Location REAL,
              PValue REAL
              )""")

    # Create a GWAS_Object tied to the empty data
    gwas = GWAS_Object(conn)

    # Test that this method raises an error
    with pytest.raises(DataNotFoundError):
        gwas.Print_QQ_Plot()

    # Because of the nature of this method, we are unsure of how to
    # test the remaining parts

    # Close the connection
    conn.close()


def test_significant_results() -> None:
    """Testing the Significant_Results method of the GWAS class."""
    # Create dummy data for testing
    data_to_insert = [
        ("Marker1", 1, 100, 0.1),
        ("Marker2", 2, 200, 0.001),
    ]

    # Create an in-memory connection with sqlite3
    conn = sqlite3.connect(":memory:")
    c = conn.cursor()

    # Replace any existing table with an empty table
    c.execute("""CREATE TABLE IF NOT EXISTS gwas (
                MarkerID TEXT,
                Chromosome TEXT,
                Location REAL,
                PValue REAL
                )""")

    # Add the dummy data
    c.executemany(
        """INSERT OR IGNORE INTO gwas 
        (MarkerID, Chromosome, Location, PValue) 
        VALUES (?, ?, ?, ?)""",
        data_to_insert,
    )

    # Commit changes
    conn.commit()

    gwas = GWAS_Object(conn)
    results = gwas.Significant_Results()

    assert results == ["Marker2"]


def test_benjamini_hochberg() -> None:
    """Testing the Benjamini_Hochberg_Procedure() function."""
    # Dummy data
    data = [
        ("A", 0.05),
        ("B", 0.01),
        ("C", 0.3),
    ]

    # Expected results after correction
    expected_results = {
        "A": False,
        "B": True,
        "C": False,
    }

    # Run the function
    actual_results = Benjamini_Hochberg_Procedure(data)

    # Assert the results
    assert actual_results == expected_results


def test_parse_gwas_output() -> None:
    """Tests that data is properly parsed."""
    # Open the connection and parse the file
    conn = sqlite3.connect(":memory:")
    Parse_GWAS_Output("test_data.csv", ",", conn)

    # Fetch the data that has been stored
    c = conn.cursor()
    c.execute("SELECT * FROM gwas")

    result = c.fetchall()
    print(result)
    assert result == [
        ("Marker", "1", 100.0, 0.001),
    ]


# Test with invalid q value
def test_benjamini_hochberg_invalid_q() -> None:
    """Tests that a ValueError is raised with invalid 'q' value."""
    # Dummy data
    data = [
        ("A", 0.1),
        ("B", 0.01),
        ("C", 0.3),
    ]

    with pytest.raises(ValueError):
        Benjamini_Hochberg_Procedure(data, q=1.2)
