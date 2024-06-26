"""Toolkit for Graphical GWAS Analysis in Python.

Created by Austin Allen, Julia Benendetti, and Jonathan Hui.

This Module Provides The Following Functions:
- Parse_GWAS_Output: Stores GWAS Output Into An SQLite Database.
    - The Input File Should Include ID, Chromasome, Location, and P-Value.

This Module Provides The Following Classes:
- GWAS_Object: An Object Containing The Output Of A GWAS Analysis.
- DataNotFoundError: An Error Called When GWAS Output Cannot Be Found.

This Module Provides The Following Attributes That A GWAS_Object Inherits:
- N_of_Significant_Tests: Provides The Number Of Sig. Tests For The Given Data.
- Alpha_Level: Determines The Significance Level With A Bonferroni Correction.
- Manhattan_Plot: Creates A Manhattan Plot For A Given GWAS Object.
- QQ_Plot: Creates A QQ Plot For A Given GWAS Object.
"""

import sqlite3
from typing import Any, Dict

import matplotlib.pyplot as plt
import numpy as np


class DataNotFoundError(Exception):
    """An Error Called When The GWAS Output Cannot Be Found."""

    pass


class GWAS_Object:
    """Represents The Output Of A Finished GWAS Analysis."""

    def __init__(self, conn: sqlite3.Connection):
        """Initializes A GWAS Output Object."""
        self.connection = conn

    @property
    def Alpha_Level(self) -> float:
        """Calculates The Adjusted Significance Level From The GWAS Data."""
        c = self.connection.cursor()
        c.execute("SELECT COUNT(*) FROM gwas WHERE PValue < ?", (0.05,))
        significant_tests_count = c.fetchone()[0]
        c.close()
        alpha: float = 0.05 / significant_tests_count
        return alpha

    def Print_Manhattan_Plot(self) -> None:
        """Creates A Manhattan Plot To Visualize The GWAS Analysis."""
        # Step 1: Retrieve Data From The SQLite Database.
        c = self.connection.cursor()
        c.execute("SELECT Chromosome, Location, PValue FROM gwas")
        data = c.fetchall()
        if not data:
            raise DataNotFoundError("No data to plot.")
        c.close()

        # Step 2: Build The Plot From Chromosome Location and P-Value.
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

        # Step 3: Initialize The Manhattan Plot.
        plt.figure(figsize=(12, 6))

        # Step 4: Stratify Each Chromosome Into Their Own Location Bins.
        label_at_bin_median = []
        current_position = 0
        for chromie, positions in chromosomes.items():
            plt.scatter(
                np.array(positions["x"]) + current_position,
                positions["y"],
                label=f"Chromosome {chromie}",
            )
            current_position += max(positions["x"])
            label_at_bin_median.append(
                current_position - max(positions["x"]) / 2
            )

        # Step 5: Draw A Horizontal Line At The Significance Level
        significancelevel: float = self.Alpha_Level
        plt.axhline(
            -np.log10(significancelevel),
            color="r",
            linestyle="--",
            label=f"-log10({significancelevel}) significance level",
        )

        # Step 6: Put Everything Together And Construct The Plot.
        plt.xlabel("Genomic Position")
        plt.ylabel("-log10(P-Value)")
        plt.title("Manhattan Plot")
        plt.xticks(
            label_at_bin_median, [f"{chromie}" for chromie in chromosomes]
        )
        plt.ylim(0)
        plt.show()

    def Print_QQ_Plot(self) -> None:
        """Creates A QQ Plot To Visualize The GWAS Analysis."""
        # Step 1: Retrieve Data From The SQLite Database.
        c = self.connection.cursor()
        c.execute("SELECT PValue FROM gwas")
        pvalues = np.array(c.fetchall(), dtype=float)
        if not pvalues.any():
            raise DataNotFoundError("No data to plot.")
        c.close()

        # Step 2: Get Theoretical Quantiles With The GWAS Output's PValues.
        n = len(pvalues)
        theoretical_quantiles = -np.log10(np.arange(1, n + 1) / n)
        plt.figure(figsize=(6, 6))
        plt.scatter(
            theoretical_quantiles,
            -np.log10(np.sort(pvalues.flatten())),
            color="blue",
            alpha=0.7,
        )

        # Step 3: Put Everything Together And Construct The Plot.
        plt.plot(
            [0, max(theoretical_quantiles)],
            [0, max(theoretical_quantiles)],
            color="red",
            linestyle="--",
        )
        plt.xlabel("Expected -log10(P-Value)")
        plt.ylabel("Observed -log10(P-Value)")
        plt.title("QQ Plot")
        plt.show()


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


# Here Is Some Sample Code To Demonstrate How This Toolkit Works!
if __name__ == "__main__":
    # Open The Connection
    connection = sqlite3.connect(":memory:")

    # Parses The Dataset
    Parse_GWAS_Output("src/demo_GWAS_Output.csv", ",", connection)
    gwasoutput = GWAS_Object(connection)

    # Find the Alpha Level, Corrected For Multiple Testing
    print(gwasoutput.Alpha_Level)

    # Print A Manhattan Plot (Warning: Pretty Colors)
    gwasoutput.Print_Manhattan_Plot()

    # Print A QQ Plot
    gwasoutput.Print_QQ_Plot()

    # Close The Connection
    connection.close()
