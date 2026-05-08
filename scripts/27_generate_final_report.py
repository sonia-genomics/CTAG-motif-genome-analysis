
# Automatically generate a final project summary
# report from all analyses.
#
# =========================================================

import os
import pandas as pd

summary_file = "output/comparative_summary.csv"

output_report = "results/CTAG_final_report.txt"

os.makedirs("results", exist_ok=True)

df = pd.read_csv(summary_file)

# BASIC STATISTICS
num_genomes = len(df)

mean_gc = round(df["GC_Content"].mean(), 2)

mean_ctag = round(df["CTAG_OE"].mean(), 3)

strong = len(
    df[df["CTAG_Category"] ==
       "Strong_Underrepresentation"]
)

moderate = len(
    df[df["CTAG_Category"] ==
       "Moderate_Underrepresentation"]
)

# GENERATE REPORT
with open(output_report, "w") as out:

    out.write("=========================================\n")
    out.write("CTAG MOTIF GENOME ANALYSIS REPORT\n")
    out.write("=========================================\n\n")

    out.write(f"Total genomes analyzed: {num_genomes}\n")
    out.write(f"Average GC content: {mean_gc}%\n")
    out.write(f"Average CTAG O/E ratio: {mean_ctag}\n\n")

    out.write("CTAG Suppression Categories\n")
    out.write("---------------------------------\n")
    out.write(f"Strong under-representation: {strong}\n")
    out.write(f"Moderate under-representation: {moderate}\n\n")

    out.write("Biological Interpretation\n")
    out.write("---------------------------------\n")

    out.write(
        "CTAG motifs show widespread under-"
        "representation across bacterial genomes. "
        "Localized enrichment regions suggest "
        "possible biological roles associated with "
        "sequence heterogeneity, genome organization, "
        "and transcription-related processes.\n\n"
    )

    out.write("Pipeline Features\n")
    out.write("---------------------------------\n")

    out.write("- Genome-wide tetranucleotide analysis\n")
    out.write("- Comparative genomics\n")
    out.write("- CTAG hotspot detection\n")
    out.write("- Statistical analysis\n")
    out.write("- PCA and clustering\n")
    out.write("- Entropy analysis\n")
    out.write("- Genome browser track generation\n")
    out.write("- Workflow reproducibility\n")

print("\n=================================================")
print("FINAL REPORT GENERATED")
print("=================================================")