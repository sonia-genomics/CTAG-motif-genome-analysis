#!/usr/bin/env python

# =========================================================
# Comparative Phylogenetic Analysis of CTAG Genomic Features
# =========================================================

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

metadata_file = "requirements/bacterial_genome.csv"

correlation_file = (
    "results/correlation_analysis/correlation_summary.csv"
)

heterogeneity_dir = "results/base_heterogeneity"

entropy_dir = "results/entropy_analysis"

output_dir = "results/phylogenetic_analysis"

os.makedirs(output_dir, exist_ok=True)

print("\nLoading metadata...")

meta = pd.read_csv(metadata_file)

print(meta.head())


print("\nLoading correlation analysis...")

corr = pd.read_csv(correlation_file)

print(corr.head())

# COLLECT ADDITIONAL FEATURES

extra_features = []

print("\nExtracting genome-level features...")

for file in os.listdir(heterogeneity_dir):

    if file.endswith("_heterogeneity.csv"):

        organism = file.replace("_heterogeneity.csv", "")

        hetero_path = os.path.join(
            heterogeneity_dir,
            file
        )

        entropy_path = os.path.join(
            entropy_dir,
            organism + "_entropy.csv"
        )

        hetero_df = pd.read_csv(hetero_path)

        # BASIC GENOMIC FEATURES
        mean_gc = hetero_df["GC_percent"].mean()

        mean_cv = hetero_df["CV"].mean()

        mean_ctag = hetero_df["CTAG_count"].mean()

        max_ctag = hetero_df["CTAG_count"].max()

        # ENTROPY FEATURES
        if os.path.exists(entropy_path):

            entropy_df = pd.read_csv(entropy_path)

            mean_entropy = entropy_df["Entropy"].mean()

            max_entropy = entropy_df["Entropy"].max()

        else:

            mean_entropy = 0
            max_entropy = 0

        extra_features.append([

            organism,

            mean_gc,
            mean_cv,

            mean_ctag,
            max_ctag,

            mean_entropy,
            max_entropy

        ])

# FEATURE DATAFRAME
feature_df = pd.DataFrame(extra_features, columns=[

    "Organism",

    "Mean_GC",
    "Mean_CV",

    "Mean_CTAG",
    "Max_CTAG",

    "Mean_Entropy",
    "Max_Entropy"

])

print("\nFeature summary:")
print(feature_df.head())

# MERGE ALL DATA
print("\nMerging datasets...")

merged = pd.merge(
    meta,
    corr,
    on="Organism"
)

merged = pd.merge(
    merged,
    feature_df,
    on="Organism"
)

print("\nMerged dataset:")
print(merged.head())


master_output = os.path.join(
    output_dir,
    "master_phylogenetic_dataset.csv"
)

merged.to_csv(master_output, index=False)

print("\nSaved:", master_output)

# SELECT FEATURES
feature_columns = [

    "Pearson_GC",
    "Pearson_CV",

    "Spearman_GC",
    "Spearman_CV",

    "Mean_GC",
    "Mean_CV",

    "Mean_CTAG",
    "Max_CTAG",

    "Mean_Entropy",
    "Max_Entropy"

]

X = merged[feature_columns]

# STANDARDIZATION
print("\nStandardizing features...")

scaler = StandardScaler()

X_scaled = scaler.fit_transform(X)

# PCA ANALYSIS
print("\nPerforming PCA...")

pca = PCA(n_components=2)

components = pca.fit_transform(X_scaled)

merged["PC1"] = components[:, 0]
merged["PC2"] = components[:, 1]

# EXPLAINED VARIANCE
explained = pca.explained_variance_ratio_

print("\nExplained variance:")
print(explained)

# SAVE PCA TABLE
pca_output = os.path.join(
    output_dir,
    "pca_coordinates.csv"
)

merged.to_csv(pca_output, index=False)

# PCA VISUALIZATION
print("\nGenerating PCA plot...")

plt.figure(figsize=(12, 10))

sns.scatterplot(

    data=merged,

    x="PC1",
    y="PC2",

    hue="family",


    s=180

)

# LABEL POINTS
for i, row in merged.iterrows():

    plt.text(

        row["PC1"],
        row["PC2"],

        row["Organism"],

        fontsize=8

    )

# TITLES
plt.title(

    "Phylogenetic CTAG Genome Architecture Analysis",

    fontsize=16,
    fontweight="bold"

)

plt.xlabel(
    f"PC1 ({explained[0]*100:.2f}% variance)"
)

plt.ylabel(
    f"PC2 ({explained[1]*100:.2f}% variance)"
)

plt.grid(alpha=0.3)

plt.tight_layout()

pca_plot = os.path.join(
    output_dir,
    "phylogenetic_pca.png"
)

plt.savefig(
    pca_plot,
    dpi=600
)

plt.close()

print("Saved:", pca_plot)

# HIERARCHICAL CLUSTERING
print("\nGenerating dendrogram...")

distance_matrix = pdist(X_scaled)

linked = linkage(
    distance_matrix,
    method='ward'
)

plt.figure(figsize=(14, 10))

dendrogram(

    linked,

    labels=merged["Organism"].tolist(),

    leaf_rotation=90,

    leaf_font_size=9

)

plt.title(

    "Hierarchical Clustering of Bacterial Genomes\nBased on CTAG Genomic Features",

    fontsize=16,
    fontweight="bold"

)

plt.ylabel("Euclidean Distance")

plt.tight_layout()

dendrogram_plot = os.path.join(
    output_dir,
    "phylogenetic_dendrogram.png"
)

plt.savefig(
    dendrogram_plot,
    dpi=600,
    bbox_inches='tight'
)

plt.close()

print("Saved:", dendrogram_plot)

# HEATMAP
print("\nGenerating heatmap...")

heatmap_data = merged.set_index("Organism")[
    feature_columns
]

plt.figure(figsize=(14, 10))

sns.heatmap(

    heatmap_data,

    cmap="coolwarm",

    annot=True,

    linewidths=0.5

)

plt.title(

    "Comparative CTAG Genomic Feature Heatmap",

    fontsize=16,
    fontweight="bold"

)

plt.tight_layout()

heatmap_plot = os.path.join(
    output_dir,
    "phylogenetic_heatmap.png"
)

plt.savefig(
    heatmap_plot,
    dpi=600,
    bbox_inches='tight'
)

plt.close()

print("Saved:", heatmap_plot)

# GROUP COMPARISON PLOTS
print("\nGenerating group comparison plots...")

comparison_features = [

    "Mean_CTAG",
    "Mean_CV",
    "Mean_Entropy",
    "Mean_GC"

]

for feature in comparison_features:

    plt.figure(figsize=(10, 6))

    sns.boxplot(

        data=merged,

        x="family",

        y=feature

    )

    sns.stripplot(

        data=merged,

        x="family",

        y=feature,

        color="black",

        alpha=0.7

    )

    plt.title(

        f"{feature} Across Bacterial Groups",

        fontsize=14,
        fontweight="bold"

    )

    plt.xticks(rotation=15)

    plt.tight_layout()

    out_path = os.path.join(

        output_dir,

        f"{feature}_group_comparison.png"

    )

    plt.savefig(
        out_path,
        dpi=600
    )

    plt.close()

    print("Saved:", out_path)

# CORRELATION MATRIX
print("\nGenerating correlation matrix...")

corr_matrix = merged[feature_columns].corr()

plt.figure(figsize=(12, 10))

sns.heatmap(

    corr_matrix,

    annot=True,

    cmap="viridis",

    linewidths=0.5

)

plt.title(

    "Correlation Matrix of CTAG Genomic Features",

    fontsize=16,
    fontweight="bold"

)

plt.tight_layout()

corr_plot = os.path.join(
    output_dir,
    "feature_correlation_matrix.png"
)

plt.savefig(
    corr_plot,
    dpi=600,
    bbox_inches='tight'
)

plt.close()

print("Saved:", corr_plot)

# FINAL SUMMARY
print("\n=================================================")
print("SUCCESS: PHYLOGENETIC CTAG ANALYSIS COMPLETED")
print("=================================================")

print("\nOutputs generated in:")
print(output_dir)

print("\nGenerated files:")

generated = os.listdir(output_dir)

for file in generated:
    print(" -", file)

print("\n=================================================")