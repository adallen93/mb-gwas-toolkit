"""GWAS Toolkit for Genomic Data Analysis in Python.

Created by Austin Allen, Julia Benendetti, and Jonathan Hui.

This Module Provides The Following Functions:
- Parse_GWAS_Output: Stores GWAS Output Into An SQLite Database.
    - The Input File Should Include ID, Chromasome, Location, and P-Value.
- Alpha_Level: Determines The Significance Level With A Bonferroni Correction.
- Manhattan_Plot: Creates A Manhattan Plot For A Given GWAS Object.
- QQ_Plot

This Module Provides The Following Classes:
- GWAS_Object: An Object Containing The Output Of A GWAS Analysis.
- MissingDataError: Called When The Genomic Data File Contains Missing Data
"""

import sqlite3
from typing import Any, Dict, Tuple


class GWAS_Object:
    """Represents A Finished GWAS Analysis."""

    def __init__(self, conn: sqlite3.Connection):
        """Initializes A GWAS Output Object."""
        self.connection = conn

    @property
    def genetic_mapping(self) -> Dict[str, Tuple[str, float]]:
        """Retrieve Chromosome and Position For Each MarkerID."""
        c = self.connection.cursor()
        c.execute("SELECT MarkerID, Chromosome, Location FROM gwas")
        results = c.fetchall()
        genemap = {result[0]: (result[1], result[2]) for result in results}
        return genemap

    @property
    def number_of_significant_tests(Self) -> Any:
        """Finds The Number Of Significant Tests."""
        # Implement A Method By Package To Calculate # Of Sig. Tests.

    def Alpha_Level(self) -> Any:
        """Calculates The Adjusted Significance Level From The GWAS Data."""
        alpha = 0.05 / self.number_of_significant_tests
        return alpha

    def Manhattan_Plot(self) -> Any:
        """Creates A Manhattan Plot To Visualize The GWAS Analysis."""
        # Implement A Method By Package To Create A Manhattan Plot.
        # Leverage The Genetic Markers Dictionary.

    def QQ_Plot(self) -> Any:
        """Creates A QQ Plot To Visualize The GWAS Analysis."""
        # Implement A Method By Package To Create A QQ Plot.
        # Leverage The Genetic Markers Dictionary.


def Parse_GWAS_Output(
    gwas_output_file: str, delimiter: str, conn: sqlite3.Connection
) -> None:
    """Parses Patient And Lab Data From Files Into An SQL Database."""
    c = conn.cursor()
    c.execute("""CREATE TABLE IF NOT EXISTS gwas (
                 MarkerID TEXT,
                 Chromosome TEXT,
                 Location REAL,
                 Phenotype TEXT,
                 FOREIGN KEY (MarkerID) REFERENCES gwas(MarkerID))""")
    with open(gwas_output_file, encoding="utf-8") as gwas_output:
        gwas_output_head = gwas_output.readline().strip().split(delimiter)
        for line in gwas_output_file:
            gwas_output_data = dict(
                zip(gwas_output_head, line.strip().split(delimiter))
            )
            c.execute(
                """INSERT OR IGNORE INTO patients 
                (MarkerID, Chromosome, Location, Phenotype)
                VALUES (?, ?, ?, ?)""",
                (
                    gwas_output_data.get("MarkerID", ""),
                    gwas_output_data.get("Chrmosome", ""),
                    gwas_output_data.get("Location", ""),
                    gwas_output_data.get("Phenotype", ""),
                ),
            )
    conn.commit()