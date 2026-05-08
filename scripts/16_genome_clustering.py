import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler
import os

INPUT_FILE = "output/final_tetra_analysis.xlsx"

OUTPUT_FIG = "results/figures/genome_dendrogram.png"

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

# SCALE
scaler = StandardScaler()

scaled = scaler.fit_transform(matrix)


# HIERARCHICAL CLUSTERING
Z = linkage(
    scaled,
    method='ward'
)


# DENDROGRAM
plt.figure(figsize=(14, 10))

dendrogram(
    Z,
    labels=labels,
    leaf_rotation=90
)

plt.title("Hierarchical Clustering of Bacterial Genomes")

plt.ylabel("Distance")

plt.tight_layout()

plt.savefig(OUTPUT_FIG, dpi=300)

print("Dendrogram saved:", OUTPUT_FIG)