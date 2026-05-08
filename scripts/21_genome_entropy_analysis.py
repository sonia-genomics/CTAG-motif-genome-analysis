import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from math import log2

input_folder = "data/processed/genomes_singleline"
output_dir = "results/entropy_analysis"

os.makedirs(output_dir, exist_ok=True)

WINDOW_SIZE = 5000
STEP_SIZE = 1000

# SHANNON ENTROPY
def shannon_entropy(sequence):

    probs = []

    for base in ["A", "T", "G", "C"]:

        p = sequence.count(base) / len(sequence)

        if p > 0:
            probs.append(p)

    entropy = -sum(p * log2(p) for p in probs)

    return entropy

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

            entropy = shannon_entropy(window)

            ctag_count = window.count("CTAG")

            gc = ((window.count("G") + window.count("C")) / len(window)) * 100

            results.append([
                start,
                end,
                entropy,
                gc,
                ctag_count
            ])

        df = pd.DataFrame(results, columns=[
            "Start",
            "End",
            "Entropy",
            "GC_percent",
            "CTAG_count"
        ])

        # SAVE CSV
        csv_path = os.path.join(
            output_dir,
            f"{organism}_entropy.csv"
        )

        df.to_csv(csv_path, index=False)

        # PLOT
        fig, ax1 = plt.subplots(figsize=(18, 6))

        ax1.plot(
            df["Start"],
            df["Entropy"],
            linewidth=1,
            label="Entropy"
        )

        ax1.set_ylabel("Shannon Entropy")

        ax2 = ax1.twinx()

        ax2.plot(
            df["Start"],
            df["CTAG_count"],
            linewidth=1,
            alpha=0.7,
            label="CTAG"
        )

        ax2.set_ylabel("CTAG Count")

        plt.title(
            f"Genome Entropy vs CTAG Density: {organism}",
            fontsize=14,
            fontweight="bold"
        )

        fig.tight_layout()

        fig_path = os.path.join(
            output_dir,
            f"{organism}_entropy_plot.png"
        )

        plt.savefig(fig_path, dpi=600)
        plt.close()

print("\nSUCCESS: Genome entropy analysis completed")
