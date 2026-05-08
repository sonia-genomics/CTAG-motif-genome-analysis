
# Generate browser-compatible BEDGRAPH tracks for:
# 1. CTAG density
# 2. GC content
# 3. Sequence entropy
#
# Lightweight version for 2 CPU / 8 GB RAM systems
# =========================================================

import os
import math
from Bio import SeqIO

genome_folder = "data/raw/genomes"

output_dir = "output/browser_tracks"
os.makedirs(output_dir, exist_ok=True)


window_size = 5000
step_size = 1000

# FUNCTIONS
def calculate_gc(seq):

    g = seq.count("G")
    c = seq.count("C")

    return ((g + c) / len(seq)) * 100 if len(seq) > 0 else 0


def calculate_entropy(seq):

    length = len(seq)

    if length == 0:
        return 0

    entropy = 0

    for base in ["A", "T", "G", "C"]:

        p = seq.count(base) / length

        if p > 0:
            entropy -= p * math.log2(p)

    return entropy


# PROCESS GENOMES
for genome_file in os.listdir(genome_folder):

    if not genome_file.endswith(".fna"):
        continue

    organism = genome_file.replace(".fna", "")

    print(f"\nProcessing: {organism}")

    genome_path = os.path.join(genome_folder, genome_file)

    record = SeqIO.read(genome_path, "fasta")

    seq = str(record.seq).upper()

    chrom = record.id


    ctag_file = open(
        os.path.join(output_dir,
        organism + "_ctag_density.bedgraph"), "w"
    )

    gc_file = open(
        os.path.join(output_dir,
        organism + "_gc_content.bedgraph"), "w"
    )

    entropy_file = open(
        os.path.join(output_dir,
        organism + "_entropy.bedgraph"), "w"
    )

    # SLIDING WINDOW ANALYSIS
    for start in range(0, len(seq) - window_size, step_size):

        end = start + window_size

        window_seq = seq[start:end]

        # CTAG DENSITY
        ctag_count = window_seq.count("CTAG")

        ctag_density = ctag_count / (window_size / 1000)

        # GC CONTENT
        gc = calculate_gc(window_seq)

        # ENTROPY
        entropy = calculate_entropy(window_seq)

        # WRITE BEDGRAPH
        ctag_file.write(
            f"{chrom}\t{start}\t{end}\t{ctag_density}\n"
        )

        gc_file.write(
            f"{chrom}\t{start}\t{end}\t{gc:.2f}\n"
        )

        entropy_file.write(
            f"{chrom}\t{start}\t{end}\t{entropy:.4f}\n"
        )

    ctag_file.close()
    gc_file.close()
    entropy_file.close()

    print("Tracks generated")

print("\n=================================================")
print("GENOME BROWSER TRACK GENERATION COMPLETED")
print("=================================================")