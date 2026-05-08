import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, kruskal
import os

INPUT_FILE = "output/final_tetra_analysis.xlsx"
OUTPUT_FILE = "results/tables/statistical_summary.csv"

os.makedirs("results/tables", exist_ok=True)

xls = pd.ExcelFile(INPUT_FILE)

results = []

for sheet in xls.sheet_names:

    df = pd.read_excel(xls, sheet_name=sheet)

    organism = sheet


    # CTAG O/E VALUES
    ctag = df[df["Tetranucleotide"] == "CTAG"]

    if ctag.empty:
        continue

    genome_ctag = ctag["Genome_OE"].values[0]
    coding_ctag = ctag["Coding_OE"].values[0]
    term_ctag = ctag["Term_OE"].values[0]


    # Mann Whitney Tests
    genome_all = df["Genome_OE"].dropna()
    coding_all = df["Coding_OE"].dropna()
    term_all = df["Term_OE"].dropna()

    try:
        stat_gc, p_gc = mannwhitneyu(genome_all, coding_all)
    except:
        stat_gc, p_gc = np.nan, np.nan

    try:
        stat_gt, p_gt = mannwhitneyu(genome_all, term_all)
    except:
        stat_gt, p_gt = np.nan, np.nan


    # DTAG KRUSKAL TEST
    motifs = ["CTAG", "GTAG", "ATAG", "TTAG"]

    values = []

    for motif in motifs:

        sub = df[df["Tetranucleotide"] == motif]

        if not sub.empty:
            values.append(sub["Genome_OE"].values)

    try:
        kruskal_stat, kruskal_p = kruskal(*values)
    except:
        kruskal_stat, kruskal_p = np.nan, np.nan

    results.append({
        "Organism": organism,

        "CTAG_Genome_OE": genome_ctag,
        "CTAG_Coding_OE": coding_ctag,
        "CTAG_Term_OE": term_ctag,

        "Genome_vs_Coding_p": p_gc,
        "Genome_vs_Termination_p": p_gt,

        "DTAG_Kruskal_p": kruskal_p
    })

summary_df = pd.DataFrame(results)

summary_df.to_csv(OUTPUT_FILE, index=False)

print("Statistical analysis completed")
print("Saved:", OUTPUT_FILE)