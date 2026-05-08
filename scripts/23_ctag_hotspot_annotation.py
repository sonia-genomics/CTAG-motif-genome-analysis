# =========================================================
# 23_ctag_hotspot_annotation.py
# =========================================================
# PURPOSE:
# Detect CTAG hotspot regions in bacterial genomes
# and annotate nearby genes/features from GenBank.
#
# OUTPUT:
# 1. CTAG hotspot coordinates
# 2. Local CTAG density
# 3. Nearby genes
# 4. Functional annotations
#
# AUTHOR: Sonia
# =========================================================

import os
import pandas as pd
from Bio import SeqIO
from collections import Counter

genome_folder = "data/raw/genomes"
annotation_folder = "data/raw/annotations"

output_dir = "output/ctag_hotspots"
os.makedirs(output_dir, exist_ok=True)

window_size = 5000
step_size = 1000

ctag_threshold = 5

# FUNCTION: FIND CTAG HOTSPOTS

def find_ctag_hotspots(sequence, window_size, step_size):

    hotspots = []

    for start in range(0, len(sequence) - window_size, step_size):

        end = start + window_size

        window_seq = sequence[start:end]

        ctag_count = window_seq.count("CTAG")

        density = ctag_count / (window_size / 1000)

        hotspots.append({
            "Start": start,
            "End": end,
            "CTAG_Count": ctag_count,
            "Density_per_kb": density
        })

    return pd.DataFrame(hotspots)

# PROCESS EACH GENOME
for genome_file in os.listdir(genome_folder):

    if not genome_file.endswith(".fna"):
        continue

    organism = genome_file.replace(".fna", "")

    print(f"\nProcessing: {organism}")

    genome_path = os.path.join(genome_folder, genome_file)
    gb_path = os.path.join(annotation_folder, organism + ".gb")

    if not os.path.exists(gb_path):
        print("GenBank file missing")
        continue


    genome_record = SeqIO.read(genome_path, "fasta")
    genome_seq = str(genome_record.seq).upper()

    # FIND HOTSPOTS
    hotspot_df = find_ctag_hotspots(
        genome_seq,
        window_size,
        step_size
    )

    # FILTER HIGH CTAG REGIONS
    hotspot_df = hotspot_df[
        hotspot_df["CTAG_Count"] >= ctag_threshold
    ].copy()

    if hotspot_df.empty:
        print("No significant hotspots found")
        continue

    # LOAD GENBANK FEATURES
    gb_record = SeqIO.read(gb_path, "genbank")

    annotations = []

    for _, row in hotspot_df.iterrows():

        hotspot_start = row["Start"]
        hotspot_end = row["End"]

        nearby_genes = []

        for feature in gb_record.features:

            if feature.type != "CDS":
                continue

            gene_start = int(feature.location.start)
            gene_end = int(feature.location.end)

            overlap = (
                gene_start <= hotspot_end and
                gene_end >= hotspot_start
            )

            if overlap:

                gene_name = feature.qualifiers.get(
                    "gene",
                    ["unknown"]
                )[0]

                product = feature.qualifiers.get(
                    "product",
                    ["unknown"]
                )[0]

                locus = feature.qualifiers.get(
                    "locus_tag",
                    ["unknown"]
                )[0]

                nearby_genes.append(
                    f"{gene_name} ({product})"
                )

        annotations.append("; ".join(nearby_genes))

    hotspot_df["Nearby_Genes"] = annotations

    output_file = os.path.join(
        output_dir,
        organism + "_ctag_hotspots.csv"
    )

    hotspot_df.to_csv(output_file, index=False)

    print("Saved:", output_file)

print("\n=================================================")
print("CTAG HOTSPOT ANNOTATION COMPLETED")
print("=================================================")