"""Demo For The GWAS Toolkit."""
# Run The Following Code To See What gwas_toolkit Does.

# Remember To First Run The Following Command In A Bash Terminal:
# pip install git+https://github.com/adallen93/mb-gwas-toolkit@main

import sqlite3

from gwas_toolkit import GWAS_Object, Parse_GWAS_Output

# Parses The Dataset
conn = sqlite3.connect(":memory:")
Parse_GWAS_Output("demo/demo_GWAS_Output.csv", ",", conn)
gwas_object = GWAS_Object(conn)

# Find the Alpha Level, Corrected For Multiple Testing
print(gwas_object.Alpha_Level)

# Print A Manhattan Plot (Warning: Pretty Colors)
gwas_object.Print_Manhattan_Plot()

# Print A QQ Plot
gwas_object.Print_QQ_Plot()
