# Python Toolkit for Graphical GWAS Representations

## Creators

* Austin Allen
* Julia Bennedetti
* Jonathan Hui


## Background
Genome wide association studies (GWAS) are observational studies that use genomes to find genetic markers for an outcome. Most commonly, GWAs are used to determine genetic markers for disease. This is done by taking genomes sequenced from subjects with and without disease, and then compared to determine genetic association. In addition to the afore mentioned, another commonly used term that might be useful is a Single nucleotide polymorphism or SNP. These are the building blocks of the DNA sequence that is being tested.  We seek to provide a toolkit to create relevant graphical representations of GWAS output, to supplement a meaningful analysis.


## For End Users:
To use this toolkit, clone this repository to your local machine running Python 3.5 or later. This toolkit requires sqlite3.  

To install our toolkit, run the following code in your terminal:

```bash
pip install git+https://github.com/adallen93/mb-gwas-toolkit@main
```

Once library has been installed, you can import it into any .py file by typing:

```python
import gwas_toolkit_v1 as gwtk
```

Remember also to import SQLite3

Before using our source code for analysis, users must be ready to provide complete results from a previously preformed GWAS. This means the data should have location (on the genome), allele, and p-value resulting from a completed GWAS analysis.  A common way to get this type of data is by using PLINK to run a GWAS analysis, and then exporting the output of that analysis as a comma-separated value (csv) file.  If the user provides exported GWAS output from PLINK, the csv file will have a column for genetic markers called MarkerID, a column for chromosome number called Chromosome, a column for location called Location, and a column for p-values called PValue.

#### Parsing GWAS Data
The package then provides a method to parse a csv file of those specifications into an SQLite database.

```python
conn = sqlite3.connect(":memory:")
gwtk.Parse_GWAS_Output("GWAS_Output.csv", ",", conn)
gwas_object = gwtk.GWAS_Object(conn)
```

#### Manhattan Plot
The following function can be used to create a Manhattan Plot from the data provided:

```python
gwas_object.Print_Manhattan_Plot()
```

Named after the NYC skyline, an ideal Manhattan Plot had chromosomes separated vertically and clearly defined peaks. The points with the highest peaks show the location of the alleles with the highest p-value, meaning they have a larger association with the observed outcome.

#### QQ-Plot
The following function will produce a QQ-Plot of the data provided:

```python
gwas_object.Print_QQ_Plot()
```
The QQ plot is a graphical representation of the deviation of the observed P values from the null hypothesis. The observed P values for each SNP are sorted from largest to smallest and plotted against expected values from a theoretical chi-squared distribution. 

Some sample data is provided in: [demo_GWAS_Output.csv](src/demo_GWAS_Output.csv).


## For Contributors: 

To run tests locally, follow these steps:

1. Clone this repository to your local machine.
2. Install the required dependencies by running: `pip install -r requirements-test.txt`
3. Run tests by running: `pytest tests/`

Make sure all tests pass before making any contributions.
