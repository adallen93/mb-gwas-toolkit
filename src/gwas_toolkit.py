"""GWAS Toolkit for Genomic Data Analysis in Python.

Created by Austin Allen, Julia Benendetti, and Jonathan Hui.

This Module Provides The Following Functions:
- Parse_GWAS_Output: Stores GWAS Output Into An SQLite Database.
    - The Input File Should Include ID, Chromasome, Location, and P-Value.
- Alpha_Level: Determines The Significance Level With A Bonferroni Correction.
- Manhattan_Plot: Creates A Manhattan Plot For A Given GWAS Object.
- QQ_Plot: Creates A QQ Plot For A Given GWAS Object.

This Module Provides The Following Classes:
- GWAS_Object: An Object Containing The Output Of A GWAS Analysis.
"""

import sqlite3
from typing import Any, Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np


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
    def number_of_significant_tests(Self) -> int:
        """Finds The Number Of Significant Tests."""
        # Implement A Method By Package To Calculate # Of Sig. Tests.
        significant_tests = 3000  # ^To Be Replaced.
        return significant_tests

    def Alpha_Level(self) -> float:
        """Calculates The Adjusted Significance Level From The GWAS Data."""
        alpha = 0.05 / self.number_of_significant_tests
        return alpha

    def Manhattan_Plot(self) -> Any:
        """Creates A Manhattan Plot To Visualize The GWAS Analysis."""
        c = self.connection.cursor()
        c.execute("SELECT Chromosome, Location, PValue FROM gwas")
        data = c.fetchall()
        if not data:
            raise ValueError("No data to plot.")
        chromosomes: Dict[Any, Any] = {}
        max_position = 0
        for chromosome, position, pvalue in data:
            chromosome = int(chromosome)
            position = int(position)
            pvalue = float(pvalue)
            if pvalue == 0:
                pvalue = 1e-300
            logpvalue = -np.log10(pvalue)
            if chromosome not in chromosomes:
                chromosomes[chromosome] = {"x": [], "y": []}
            chromosomes[chromosome]["x"].append(position)
            chromosomes[chromosome]["y"].append(logpvalue)
            max_position = max(max_position, position)
        plt.figure(figsize=(12, 6))
        label_at_midpoints = []
        current_position = 0
        for chromie, positions in chromosomes.items():
            plt.scatter(
                np.array(positions["x"]) + current_position,
                positions["y"],
                label=f"Chromosome {chromie}",
            )
            current_position += max(positions["x"])
            label_at_midpoints.append(
                current_position - max(positions["x"]) / 2
            )
        plt.xlabel("Genomic Position")
        plt.ylabel("-log10(P-Value)")
        plt.title("Manhattan Plot")
        plt.xticks(
            label_at_midpoints, [f"{chromie}" for chromie in chromosomes]
        )
        plt.ylim(0)
        plt.show()

    def QQ_Plot(self) -> Any:
        """Creates A QQ Plot To Visualize The GWAS Analysis."""
        # Implement A Method By Package To Create A QQ Plot.
        # Leverage The Genetic Markers Dictionary and Get P-Values.
        # Then Plot P-Values Across Chromosome Bins By Location.


def Parse_GWAS_Output(
    gwas_output_file: str, delimiter: str, conn: sqlite3.Connection
) -> None:
    """Parses GWAS Output Files From Files Into An SQL Database."""
    c = conn.cursor()
    c.execute("""CREATE TABLE IF NOT EXISTS gwas (
                 MarkerID TEXT,
                 Chromosome TEXT,
                 Location REAL,
                 PValue TEXT,
                 FOREIGN KEY (MarkerID) REFERENCES gwas(MarkerID))""")
    with open(gwas_output_file, encoding="utf-8") as gwas_output:
        gwas_output_head = gwas_output.readline().strip().split(delimiter)
        for line in gwas_output:
            gwas_output_data = dict(
                (key.strip('"'), value.strip('"'))
                for key, value in zip(
                    gwas_output_head, line.strip().split(delimiter)
                )
            )
            c.execute(
                """INSERT or IGNORE INTO gwas
                (MarkerID, Chromosome, Location, PValue)
                VALUES (?, ?, ?, ?)""",
                (
                    gwas_output_data.get("MarkerID", ""),
                    gwas_output_data.get("Chromosome", ""),
                    gwas_output_data.get("Location", ""),
                    gwas_output_data.get("PValue", ""),
                ),
            )
    conn.commit()


conn = sqlite3.connect(":memory:")
Parse_GWAS_Output("GWAS_Output.csv", ",", conn)
gwas_obj = GWAS_Object(conn)
gwas_obj.Manhattan_Plot()
