import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from scipy.stats import poisson

input_folder = "data/processed/genomes_singleline"
output_dir = "results/ctag_local_clustering"

os.makedirs(output_dir, exist_ok=True)

# PARAMETERS
MOTIF = "CTAG"
WINDOW_SIZE = 5000
STEP_SIZE = 1000
ZOOM_THRESHOLD = 3

# FUNCTIONS
def sliding_window_ctag(seq, motif, window_size, step_size):

    positions = []
    observed_counts = []
    expected_counts = []
    oe_ratios = []

    genome_length = len(seq)

    total_ctag = seq.count(motif)

    # Poisson expectation
    expected_lambda = (total_ctag / genome_length) * window_size

    for start in range(0, genome_length - window_size, step_size):

        end = start + window_size

        window_seq = seq[start:end]

        observed = window_seq.count(motif)

        expected = expected_lambda

        oe = observed / expected if expected > 0 else 0

        positions.append(start)
        observed_counts.append(observed)
        expected_counts.append(expected)
        oe_ratios.append(oe)

    return pd.DataFrame({
        "Position": positions,
        "Observed": observed_counts,
        "Expected": expected_counts,
        "OE": oe_ratios
    })


# PROCESS EACH GENOME
for genome_file in os.listdir(input_folder):

    if genome_file.endswith(".fna"):

        organism = genome_file.replace(".fna", "")

        print(f"Processing: {organism}")

        genome_path = os.path.join(input_folder, genome_file)

        record = SeqIO.read(genome_path, "fasta")

        sequence = str(record.seq).upper()

        # SLIDING WINDOW ANALYSIS
        df = sliding_window_ctag(
            sequence,
            MOTIF,
            WINDOW_SIZE,
            STEP_SIZE
        )

        # IDENTIFY LOCAL CLUSTERS
        cluster_df = df[df["OE"] >= ZOOM_THRESHOLD]

        # TABLE
        table_path = os.path.join(
            output_dir,
            f"{organism}_local_clusters.csv"
        )

        cluster_df.to_csv(table_path, index=False)

        # MAIN FIGURE
        fig = plt.figure(figsize=(18, 10))

        gs = fig.add_gridspec(
            3,
            3,
            height_ratios=[2, 1, 1]
        )

        # TOP PANEL
        ax_main = fig.add_subplot(gs[0, :])

        ax_main.plot(
            df["Position"],
            df["Observed"],
            color="deepskyblue",
            linewidth=1,
            label="Observed CTAG"
        )

        ax_main.plot(
            df["Position"],
            df["Expected"],
            color="red",
            linewidth=1.5,
            label="Expected CTAG"
        )

        # Highlight cluster regions
        for _, row in cluster_df.iterrows():

            ax_main.axvspan(
                row["Position"],
                row["Position"] + WINDOW_SIZE,
                color="green",
                alpha=0.2
            )

        ax_main.set_title(
            f"Local Over-Representation of CTAG in {organism}",
            fontsize=16,
            fontweight="bold"
        )

        ax_main.set_ylabel("CTAG Count")

        ax_main.legend()

        ax_main.grid(alpha=0.3)

        # ZOOM PANELS
        max_zoom = min(4, len(cluster_df))

        for i in range(max_zoom):

            row = cluster_df.iloc[i]

            start = int(row["Position"])
            end = start + WINDOW_SIZE

            zoom_df = df[
                (df["Position"] >= start - 10000) &
                (df["Position"] <= end + 10000)
            ]

            ax_zoom = fig.add_subplot(gs[1 + i // 2, i % 2])

            ax_zoom.plot(
                zoom_df["Position"],
                zoom_df["Observed"],
                color="deepskyblue",
                linewidth=1
            )

            ax_zoom.plot(
                zoom_df["Position"],
                zoom_df["Expected"],
                color="red",
                linewidth=1.5
            )

            ax_zoom.axvspan(
                start,
                end,
                color="green",
                alpha=0.2
            )

            ax_zoom.set_title(
                f"Cluster {i+1}",
                fontsize=10
            )

            ax_zoom.grid(alpha=0.3)

        plt.tight_layout()

        fig_path = os.path.join(
            output_dir,
            f"{organism}_ctag_cluster_analysis.png"
        )

        plt.savefig(
            fig_path,
            dpi=600,
            bbox_inches="tight"
        )

        plt.close()

        print(f"Saved: {fig_path}")

print("\nSUCCESS: CTAG local clustering analysis completed")
