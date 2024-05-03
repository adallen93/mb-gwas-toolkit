"""Demo For The GWAS Toolkit."""
# Run The Following Code To See What gwas_toolkit Does.

# Remember To First Run The Following Command In A Bash Terminal:
# pip install git+https://github.com/adallen93/mb-gwas-toolkit@main

import sqlite3

from gwas_toolkit import GWAS_Object, Parse_GWAS_Output

# Open The Connection
connection = sqlite3.connect(":memory:")

# Parses The Dataset
Parse_GWAS_Output("demo/demo_GWAS_Output.csv", ",", connection)
gwasoutput = GWAS_Object(connection)

# Find the Alpha Level, Corrected For Multiple Testing
print(gwasoutput.Alpha_Level)

# Print A Manhattan Plot (Warning: Pretty Colors)
gwasoutput.Print_Manhattan_Plot()

# Print A QQ Plot
gwasoutput.Print_QQ_Plot()

# Close The Connection
connection.close()
