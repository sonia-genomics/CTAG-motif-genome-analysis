import os
import pandas as pd
from Bio import SeqIO


genes_folder = "output/genes"
output_path = "output/last4_nt.xlsx"


with pd.ExcelWriter(output_path, engine="openpyxl") as writer:

    for file in os.listdir(genes_folder):

        if file.endswith(".fasta") or file.endswith(".fa"):

            file_path = os.path.join(genes_folder, file)
            accession = file.replace(".fasta","").replace(".fa","")

            print("Processing:", accession)

            results = []

            for record in SeqIO.parse(file_path, "fasta"):

                header = record.description
                seq = str(record.seq)

                parts = header.split()
                gene_name = parts[0]
                coords = parts[1] if len(parts) > 1 else "NA"
                strand = parts[2].strip("()") if len(parts) > 2 else "NA"


                last3 = seq[-3:] if len(seq) >= 3 else seq
                last4 = seq[-4:] if len(seq) >= 4 else seq

                row = {
                    "Gene": gene_name,
                    "Coordinates": coords,
                    "Strand": strand,
                    "Gene_Length": len(seq),
                    "Last_3_nt": last3,
                    "Last_4_nt": last4
                }

                results.append(row)

            df = pd.DataFrame(results)

            sheet_name = accession[:31]  # Excel limit
            df.to_excel(writer, sheet_name=sheet_name, index=False)

            print("Saved sheet:", sheet_name)

print("\n DONE! File saved at:", output_path)