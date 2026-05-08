
# Create a compact comparative genomics summary
# across all bacterial genomes.
#
# Lightweight version
# =========================================================

import os
import pandas as pd


qc_file = "output/genome_qc_results.csv"
tetra_file = "output/CTAG_permutation_final.xlsx"

output_file = "output/comparative_summary.csv"

qc_df = pd.read_csv(qc_file)

genome_df = pd.read_excel(
    tetra_file,
    sheet_name="genome_permutation"
)

# EXTRACT CTAG VALUES
summary = genome_df[[
    "BACTERIAL SPECIES",
    "CTAG",
    "GTAC",
    "CATG"
]].copy()

summary.rename(columns={
    "BACTERIAL SPECIES": "Genome",
    "CTAG": "CTAG_OE",
    "GTAC": "GTAC_OE",
    "CATG": "CATG_OE"
}, inplace=True)

# MERGE QC + CTAG DATA
final_df = pd.merge(
    qc_df,
    summary,
    on="Genome",
    how="left"
)

# CLASSIFY CTAG SUPPRESSION
def classify_ctag(x):

    if x < 0.5:
        return "Strong_Underrepresentation"

    elif x < 0.8:
        return "Moderate_Underrepresentation"

    else:
        return "Weak_Underrepresentation"

final_df["CTAG_Category"] = final_df[
    "CTAG_OE"
].apply(classify_ctag)


final_df.to_csv(output_file, index=False)

print("\n=================================================")
print("COMPARATIVE SUMMARY GENERATED")
print("=================================================")