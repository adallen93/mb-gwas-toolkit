import pandas as pd

# Read the original CSV file
df = pd.read_csv("GWAS_Output.csv")

# Set values greater than 1 in the "PValue" column to 1
df["PValue"] = df["PValue"].apply(lambda x: 1 if x > 1 else x)

# Write the modified DataFrame to a new CSV file
df.to_csv("GWAS_Output_For_Tests.csv", index=False)
