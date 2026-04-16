import os
from Bio import SeqIO

input_folder = "output/genes"
output_folder = "output/coding_fasta"

os.makedirs(output_folder, exist_ok=True)

for file in os.listdir(input_folder):
    if file.endswith(".fasta") or file.endswith(".fa"):
        
        input_path = os.path.join(input_folder, file)
        output_path = os.path.join(output_folder, file.replace(".fasta", "_combined.fasta").replace(".fa", "_combined.fasta"))

        combined_seq = ""

        for record in SeqIO.parse(input_path, "fasta"):
            combined_seq += str(record.seq)

        with open(output_path, "w") as out:
            out.write(f">{file}_combined\n")
            out.write(combined_seq + "\n")

        print(f"Done: {file}")