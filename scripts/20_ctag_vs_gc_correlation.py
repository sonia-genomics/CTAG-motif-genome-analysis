import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression
import numpy as np


input_dir = "results/base_heterogeneity"
output_dir = "results/correlation_analysis"

os.makedirs(output_dir, exist_ok=True)

summary_results = []

# ANALYSIS
for file in os.listdir(input_dir):

    if file.endswith("_heterogeneity.csv"):

        organism = file.replace("_heterogeneity.csv", "")

        print(f"Processing: {organism}")

        path = os.path.join(input_dir, file)

        df = pd.read_csv(path)

        # CORRELATIONS
        pearson_gc, p_gc = pearsonr(
            df["CTAG_count"],
            df["GC_percent"]
        )

        pearson_cv, p_cv = pearsonr(
            df["CTAG_count"],
            df["CV"]
        )

        spearman_gc, _ = spearmanr(
            df["CTAG_count"],
            df["GC_percent"]
        )

        spearman_cv, _ = spearmanr(
            df["CTAG_count"],
            df["CV"]
        )

        summary_results.append([
            organism,
            pearson_gc,
            p_gc,
            pearson_cv,
            p_cv,
            spearman_gc,
            spearman_cv
        ])

        # REGRESSION PLOT
        X = df[["GC_percent"]]
        y = df["CTAG_count"]

        model = LinearRegression()
        model.fit(X, y)

        pred = model.predict(X)

        plt.figure(figsize=(8, 6))

        plt.scatter(
            df["GC_percent"],
            df["CTAG_count"],
            alpha=0.5
        )

        plt.plot(
            df["GC_percent"],
            pred,
            linewidth=2
        )

        plt.xlabel("GC %")
        plt.ylabel("CTAG Count")

        plt.title(
            f"CTAG vs GC Correlation: {organism}"
        )

        plot_path = os.path.join(
            output_dir,
            f"{organism}_gc_correlation.png"
        )

        plt.savefig(plot_path, dpi=600)
        plt.close()

# SUMMARY TABLE
summary_df = pd.DataFrame(summary_results, columns=[
    "Organism",
    "Pearson_GC",
    "Pvalue_GC",
    "Pearson_CV",
    "Pvalue_CV",
    "Spearman_GC",
    "Spearman_CV"
])

summary_df.to_csv(
    os.path.join(output_dir, "correlation_summary.csv"),
    index=False
)

print("\nSUCCESS: Correlation analysis completed")
