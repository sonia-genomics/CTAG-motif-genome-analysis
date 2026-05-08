import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

INPUT_FILE = "output/final_tetra_analysis.xlsx"

OUTPUT_FIG = "results/figures/tetra_heatmap.png"

os.makedirs("results/figures", exist_ok=True)

xls = pd.ExcelFile(INPUT_FILE)

heatmap_data = []

for sheet in xls.sheet_names:

    df = pd.read_excel(xls, sheet_name=sheet)

    row = {"Organism": sheet}

    for _, r in df.iterrows():

        motif = r["Tetranucleotide"]
        oe = r["Genome_OE"]

        row[motif] = oe

    heatmap_data.append(row)

heatmap_df = pd.DataFrame(heatmap_data)

heatmap_df.set_index("Organism", inplace=True)

# CLUSTER HEATMAP
sns.clustermap(
    heatmap_df,
    cmap="viridis",
    metric="euclidean",
    method="average",
    figsize=(18, 12)
)

plt.savefig(OUTPUT_FIG, dpi=300)

print("Heatmap saved:", OUTPUT_FIG)