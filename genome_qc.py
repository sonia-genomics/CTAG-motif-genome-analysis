import os
from Bio import SeqIO
import pandas as pd

genome_folder = "output/genomes"

results = []

for file in os.listdir(genome_folder):

    if file.endswith(".fna"):

        path = os.path.join(genome_folder, file)

        for record in SeqIO.parse(path, "fasta"):

            seq = str(record.seq).upper()

            length = len(seq)

            A = seq.count("A")
            T = seq.count("T")
            G = seq.count("G")
            C = seq.count("C")
            N = seq.count("N")

            gc_content = ((G + C) / length) * 100

            results.append([
                file,
                length,
                gc_content,
                N
            ])

df = pd.DataFrame(results, columns=[
    "Genome",
    "Genome_Length",
    "GC_Content",
    "Ambiguous_N"
])

df.to_csv("output/genome_qc_results.csv", index=False)

print("Genome QC completed")
