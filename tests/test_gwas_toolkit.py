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

import io
import sqlite3
import unittest
from unittest.mock import patch

import pytest
from gwas_toolkit import (
    Benjamini_Hochberg_Procedure,
    DataNotFoundError,
    GWAS_Object,
    Parse_GWAS_Output,
)


class TestGWASAnalysis(unittest.TestCase):
    """Creating class."""

    def setUp(self) -> None:
        """Create an in-memory SQLite database and parse test data."""
        self.conn = sqlite3.connect(":memory:")
        # Prepare test data
        test_data = """MarkerID,Chromosome,Location,PValue
        Marker1,1,100,0.001
        Marker2,2,200,0.005
        Marker3,3,300,0.0001
        """
        # Call the function to parse GWAS output
        Parse_GWAS_Output(test_data, ",", self.conn)

    def test_Parse_GWAS_Output(self) -> None:
        """Testing Parsing Function."""
        # Check if the data is correctly parsed and stored in the database
        c = self.conn.cursor()
        c.execute("SELECT * FROM gwas")
        data = c.fetchall()
        self.assertEqual(len(data), 3)
        self.assertEqual(data[0], ("Marker1", "1", 100.0, "0.001"))
        self.assertEqual(data[1], ("Marker2", "2", 200.0, "0.005"))
        self.assertEqual(data[2], ("Marker3", "3", 300.0, "0.0001"))

    def test_GWAS_Object_Manhattan_Plot(self) -> None:
        """Testing Manhattan Plot Function."""
        # Prepare test data
        test_data = [
            ("1", "100", "0.001"),
            ("2", "200", "0.005"),
            ("3", "300", "0.0001"),
        ]
        c = self.conn.cursor()
        c.executemany(
            "INSERT INTO gwas (MarkerID, Chromosome, Location, PValue) VALUES (?, ?, ?, ?)",
            test_data,
        )
        self.conn.commit()

        # Instantiate GWAS_Object
        gwas_obj = GWAS_Object(self.conn)

        # Patch plt.show() to prevent plotting and capture stdout
        with patch("matplotlib.pyplot.show"), patch(
            "sys.stdout", new_callable=io.StringIO
        ) as fake_stdout:
            # Call Manhattan_Plot method
            result = gwas_obj.Manhattan_Plot

            # Check if Manhattan_Plot returns the expected message
            self.assertIn("Manhattan Plot", fake_stdout.getvalue())
            self.assertIn(
                "A Separate Window Displaying The Manhattan Plot Was Opened.",
                result,
            )

    def test_GWAS_Object_QQ_Plot(self) -> None:
        """Testing QQ-Plot FUnction."""
        # Prepare test data
        test_data = [("1", "0.001"), ("2", "0.005"), ("3", "0.0001")]
        c = self.conn.cursor()
        c.executemany(
            "INSERT INTO gwas (MarkerID, PValue) VALUES (?, ?)", test_data
        )
        self.conn.commit()

        # Instantiate GWAS_Object
        gwas_obj = GWAS_Object(self.conn)

        # Patch plt.show() to prevent plotting and capture stdout
        with patch("matplotlib.pyplot.show"), patch(
            "sys.stdout", new_callable=io.StringIO
        ) as fake_stdout:
            # Call QQ_Plot method
            result = gwas_obj.QQ_Plot

            # Check if QQ_Plot returns the expected message
            self.assertIn("QQ Plot", fake_stdout.getvalue())
            self.assertIn(
                "A Separate Window Displaying he QQ Plot Was Opened.", result
            )

    def test_GWAS_Object_Significant_Results(self) -> None:
        """Testing Significance Function."""
        # Prepare test data
        test_data = [
            ("Marker1", "0.001"),
            ("Marker2", "0.0005"),
            ("Marker3", "0.01"),
        ]
        c = self.conn.cursor()
        c.executemany(
            "INSERT INTO gwas (MarkerID, PValue) VALUES (?, ?)", test_data
        )
        self.conn.commit()

        # Instantiate GWAS_Object
        gwas_obj = GWAS_Object(self.conn)

        # Call Significant_Results method
        significant_markers = gwas_obj.Significant_Results()

        # Check if the method returns the correct significant markers
        self.assertEqual(significant_markers, ["Marker1", "Marker2"])

    def test_GWAS_Object_No_Significant_Genes_Error(self) -> None:
        """Testing No Significant Genes."""
        # Instantiate GWAS_Object with an empty database
        gwas_obj = GWAS_Object(self.conn)

        # Call Significant_Results method
        with self.assertRaises(NoSignificantGenesError):
            gwas_obj.Significant_Results()

    def test_GWAS_Object_Data_Not_Found_Error(self) -> None:
        """Testing Object Data Not Found."""
        # Instantiate GWAS_Object with an empty database
        gwas_obj = GWAS_Object(self.conn)

        # Call Manhattan_Plot method
        with self.assertRaises(DataNotFoundError):
            gwas_obj.Print_Manhattan_Plot()

        # Call QQ_Plot method
        with self.assertRaises(DataNotFoundError):
            gwas_obj.Print_QQ_Plot()


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


# Test with invalid data types
def test_benjamini_hochberg_invalid_data() -> None:
    """Tests that a TypeError is raised with invalid types in 'data'."""
    with pytest.raises(TypeError):
        Benjamini_Hochberg_Procedure(
            [1, 2, 3]
        )  # List of numbers instead of tuples


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


if __name__ == "__main__":
    unittest.main()
