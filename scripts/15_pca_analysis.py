import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os

INPUT_FILE = "output/final_tetra_analysis.xlsx"

OUTPUT_FIG = "results/figures/pca_clustering.png"

os.makedirs("results/figures", exist_ok=True)

xls = pd.ExcelFile(INPUT_FILE)

data = []

labels = []

for sheet in xls.sheet_names:

    df = pd.read_excel(xls, sheet_name=sheet)

    vector = []

    for _, row in df.iterrows():

        vector.append(row["Genome_OE"])

    data.append(vector)

    labels.append(sheet)


matrix = pd.DataFrame(data)


# STANDARDIZE
scaler = StandardScaler()

scaled = scaler.fit_transform(matrix)


# PCA
pca = PCA(n_components=2)

coords = pca.fit_transform(scaled)


# PLOT
plt.figure(figsize=(10, 8))

plt.scatter(coords[:, 0], coords[:, 1])

for i, label in enumerate(labels):

    plt.text(
        coords[i, 0],
        coords[i, 1],
        label,
        fontsize=8
    )

plt.xlabel("PC1")
plt.ylabel("PC2")

plt.title("PCA of Tetranucleotide O/E Profiles")

plt.tight_layout()

plt.savefig(OUTPUT_FIG, dpi=300)

print("PCA plot saved:", OUTPUT_FIG)