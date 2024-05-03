"""Demo For The GWAS Toolkit."""

import sqlite3

from gwas_toolkit import GWAS_Object, Parse_GWAS_Output

# Parses The Dataset
conn = sqlite3.connect(":memory:")
Parse_GWAS_Output("demo/demo_GWAS_Output.csv", ",", conn)
gwas_object = GWAS_Object(conn)

# Print A Manhattan Plot (Warning: Pretty Colors)
gwas_object.Print_Manhattan_Plot()

# Print A QQ Plot
gwas_object.Print_QQ_Plot()

# Identifying Significant Alleles
gwas_object.Significant_Results()
