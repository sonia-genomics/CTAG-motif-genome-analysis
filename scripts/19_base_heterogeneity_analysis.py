import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

input_folder = "data/processed/genomes_singleline"
cluster_folder = "results/ctag_local_clustering"
output_dir = "results/base_heterogeneity"

os.makedirs(output_dir, exist_ok=True)

# PARAMETERS
WINDOW_SIZE = 5000
STEP_SIZE = 1000

# FUNCTION
def calculate_cv(a, t, g, c):

    arr = np.array([a, t, g, c])

    mean = np.mean(arr)
    std = np.std(arr)

    if mean == 0:
        return 0

    return std / mean


# ANALYSIS
for genome_file in os.listdir(input_folder):

    if genome_file.endswith(".fna"):

        organism = genome_file.replace(".fna", "")

        print(f"Processing: {organism}")

        genome_path = os.path.join(input_folder, genome_file)

        record = SeqIO.read(genome_path, "fasta")

        seq = str(record.seq).upper()

        results = []

        for start in range(0, len(seq) - WINDOW_SIZE, STEP_SIZE):

            end = start + WINDOW_SIZE

            window = seq[start:end]

            A = window.count("A")
            T = window.count("T")
            G = window.count("G")
            C = window.count("C")

            gc = ((G + C) / len(window)) * 100

            ctag_count = window.count("CTAG")

            cv = calculate_cv(A, T, G, C)

            results.append([
                start,
                end,
                A,
                T,
                G,
                C,
                gc,
                ctag_count,
                cv
            ])

        df = pd.DataFrame(results, columns=[
            "Start",
            "End",
            "A",
            "T",
            "G",
            "C",
            "GC_percent",
            "CTAG_count",
            "CV"
        ])

        # SAVE CSV
        csv_path = os.path.join(
            output_dir,
            f"{organism}_heterogeneity.csv"
        )

        df.to_csv(csv_path, index=False)

        # VISUALIZATION
        fig, ax1 = plt.subplots(figsize=(18, 6))

        ax1.plot(
            df["Start"],
            df["CV"],
            color="darkorange",
            linewidth=1,
            label="Base Heterogeneity (CV)"
        )

        ax1.set_ylabel("CV")
        ax1.set_xlabel("Genome Position")

        ax2 = ax1.twinx()

        ax2.plot(
            df["Start"],
            df["CTAG_count"],
            color="deepskyblue",
            linewidth=1,
            label="CTAG Count"
        )

        ax2.set_ylabel("CTAG Count")

        plt.title(
            f"CTAG vs Base Heterogeneity: {organism}",
            fontsize=14,
            fontweight="bold"
        )

        fig.tight_layout()

        fig_path = os.path.join(
            output_dir,
            f"{organism}_heterogeneity_plot.png"
        )

        plt.savefig(fig_path, dpi=600)

        plt.close()

print("\nSUCCESS: Base heterogeneity analysis completed")

