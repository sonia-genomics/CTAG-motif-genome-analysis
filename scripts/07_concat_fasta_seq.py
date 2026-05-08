import os
from Bio import SeqIO

input_folder = "data/raw/genomes"
output_folder = "data/processed/genomes_singleline"

os.makedirs(output_folder, exist_ok=True)

for file in os.listdir(input_folder):

    if file.endswith(".fna") or file.endswith(".fasta") or file.endswith(".fa"):

        input_path = os.path.join(input_folder, file)
        output_path = os.path.join(output_folder, file)

        print("Processing:", file)

        with open(output_path, "w") as out:

            for record in SeqIO.parse(input_path, "fasta"):

                out.write(f">{record.description}\n")
                out.write(str(record.seq) + "\n")

        print("Saved:", output_path)

print("\nDONE! All genomes converted to single-line FASTA")