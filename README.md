# GWAS Toolkit

## Creators

* Julia Bennedetti
* Jonathan Hui
* Austin Allen

## Background
Genome wide association studies (GWAS) are observational studies that use genomes to find genetic markers for an outcome. Most commonly, GWAs are used to determine genetic markers for disease. This is done by taking genomes sequenced from subjects with and without disease, and then compared to determine genetic association. In addition to the afore mentioned, another commonly used term that might be useful is a Single nucleotide polymorphism or SNP. These are the building blocks of the DNA sequence that is being tested.

## Using Toolkit

**Installation**  
To install our toolkit, runt he following code in your terminal:

```bash
pip install git+https://github.com/adallen93/mb-gwas-toolkit@main
```

**Usage**  
Once library has been installed, you can import it into any .py file as follows:

```python
import gwas_toolkit as gwtk
```
Before using our source code for analysis, users must be ready to provide complete results from a previously preformed GWAS. This means the data should have location (on the genome), allele, and p-value resulting from the GWAS. The following functions can be used to suport findings of the GWAS or a sensetivity analysis. 

**Parsing**  
The following function can be used to Parse Data that is a Delimiter Separated Values table with a column for genetic markers called MarkerID, a column for chromosome number called Chromosome, a column for location called Location, and a column for p-values called PValue. These are consistent with PLINK GWAS Output knitted into csv/tsv files.

```python
def Parse_GWAS_Output(
    gwas_output_file: str, delimiter: str, conn: sqlite3.Connection
) -> None:
    """Parses GWAS Output Files From Files Into An SQL Database."""
    # Create a cursor to interact with the database
    c = conn.cursor()

    # Create the 'gwas' table if it doesn't already exist
    c.execute("""CREATE TABLE IF NOT EXISTS gwas (
                 MarkerID TEXT,
                 Chromosome TEXT,
                 Location REAL,
                 PValue TEXT,
                 FOREIGN KEY (MarkerID) REFERENCES gwas(MarkerID))""")

    # Open the GWAS output file
    with open(gwas_output_file, encoding="utf-8") as gwas_output:
        # Read the header line and split it into column names
        gwas_output_head = gwas_output.readline().strip().split(delimiter)

        # Iterate over each line in the GWAS output file
        for line in gwas_output:
            # Split the line and create a dictionary with column names as keys
            gwas_output_data = dict(
                (key.strip('"'), value.strip('"'))
                for key, value in zip(
                    gwas_output_head, line.strip().split(delimiter)
                )
            )

            # Insert data into the 'gwas' table, ignoring if MarkerID exists
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

    # Commit changes to the database
    conn.commit()
```

**Manhattan Plot**  
The following function can be used to create a Manhattan Plot from the data provided:

```python
 def Print_Manhattan_Plot(self) -> None:
        """Creates A Manhattan Plot To Visualize The GWAS Analysis."""
        # Step 1: Retrieve Data From The SQLite Database.
        c = self.connection.cursor()
        c.execute("SELECT Chromosome, Location, PValue FROM gwas")
        data = c.fetchall()
        if not data:
            raise DataNotFoundError("Unable to print plot.")

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

        # Print the plot and a message to the user stating the plot has printed
        plt.show()
        print("A separate window displaying the Manhattan plot was opened.")
```

Named after the NYC skyline, an ideal Manhattan Plot had chromosomes separated vertically and clearly defined peaks. The points with the highest peaks show the location of the alleles with the highest p-value, meaning they have a larger association with the observed outcome.

**QQ-Plot**  
The following function will produce a QQ-Plot of the data provided:

```python
def Print_QQ_Plot(self) -> None:
        """Creates A QQ Plot To Visualize The GWAS Analysis."""
        # Step 1: Retrieve Data From The SQLite Database.
        c = self.connection.cursor()
        c.execute("SELECT PValue FROM gwas")
        pvalues = np.array(c.fetchall(), dtype=float)
        if not pvalues.any():
            raise DataNotFoundError("Unable to print plot.")

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

        # Print the plot and a message to the user stating the plot has printed
        plt.show()
        print("A separate window displaying the Manhattan plot was opened.")

```
The QQ plot is a graphical representation of the deviation of the observed P values from the null hypothesis. The observed P values for each SNP are sorted from largest to smallest and plotted against expected values from a theoretical chi-squared distribution. 

**Significance**  
The following function was created to help users determine significant alleles from the study. GWAS can be tested using millions of Single nucleotide polymorphisms (SNPs). 

``` python
def Significant_Results(self) -> List[str]:
        """Identifies The Names Of Gene Markers With Significant P-Values."""
        # Use connection cursor to retrieve p-values and marker ID
        c = self.connection.cursor()
        c.execute("SELECT MarkerID, PValue FROM gwas")
        data = c.fetchall()

        # Determine significant tests based on the Benjamani-Hochberg procedure
        test_results = Benjamini_Hochberg_Procedure(data, self.Alpha_Level)
        significant_gene_markers: list[str] = []

        # Warn the user if no significant gene markers were found
        if not test_results:
            warning("No significant gene markers were found.")

        else:
            # Add significant results to the final
            for id, result in test_results.items():
                if result:
                    significant_gene_markers.append(id)

        return significant_gene_markers
```
**Correction**  
One of the most important parts of GWAS is multiple testing correction. Each subject provides a genome containing millions of individual observations creating correlated data. The fllowing function has been created to run a multiple testing correction using the Benjamini Hochberg Method.

```python
def Benjamini_Hochberg_Procedure(
    data: list[tuple[str, float]], q: float = 0.05
) -> dict[str, bool]:
    """Multiple testing correction which controls FDR."""
    test_results: dict[str, bool] = {}

    # Type check for `data`
    if not all(
        isinstance(item, tuple)
        and len(item) == 2
        and isinstance(item[0], str)
        and isinstance(item[1], float)
        for item in data
    ):
        raise TypeError(
            "data must be a list of tuples where each tuple is (str, float)"
        )

    # Ensure `q` is between 0 and 1
    if not 0 < q < 1:
        raise ValueError("q must be between 0 and 1")

    # p-values must be sorted for the Benjamanai-Hochberg method
    data.sort(key=lambda x: x[1])

    # Conduct the corrected test for each p-value
    for i, results in enumerate(data):
        # Extract id and p-value from results for readability
        id = results[0]
        p_value = results[1]
        index = i + 1  # per procedure requirements

        # Add the corrected test results to `test_results` under `id`
        # Note that `True` indicates there is a result, meaning H0 was rejected
        test_results[id] = p_value < (q * index) / len(data)

    return test_results

```