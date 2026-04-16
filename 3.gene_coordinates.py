import os
import pandas as pd
from Bio import SeqIO

# -------- INPUT --------
genbank_folder = "output/annotations"
output_path = "output/gene_coordinates_per_genome.xlsx"

# ---------- WRITE MULTIPLE SHEETS ----------
with pd.ExcelWriter(output_path, engine="openpyxl") as writer:

    for gb_file in os.listdir(genbank_folder):
        if gb_file.endswith(".gb"):

            accession = gb_file.replace(".gb", "")
            gb_path = os.path.join(genbank_folder, gb_file)

            print("Processing:", accession)

            genes_list = []

            for record in SeqIO.parse(gb_path, "genbank"):
                for feature in record.features:

                    if feature.type == "CDS":

                        # Gene name handling
                        gene_name = feature.qualifiers.get("gene", [""])[0]
                        if not gene_name:
                            gene_name = feature.qualifiers.get("locus_tag", ["unknown"])[0]

                        # Coordinates
                        start = int(feature.location.start) + 1
                        end = int(feature.location.end)

                        # Extra biological info
                        strand = "+" if feature.location.strand == 1 else "-"
                        length = end - start + 1

                        # Sequence + CTAG count
                        sequence = feature.extract(record.seq)
                        ctag_count = sequence.count("CTAG")

                        genes_list.append({
                            "gene": gene_name,
                            "start": start,
                            "end": end,
                            "strand": strand,
                            "length": length,
                        })
            df = pd.DataFrame(genes_list)
            sheet_name = accession[:31]
            df.to_excel(writer, sheet_name=sheet_name, index=False)

print("Multi-sheet Excel file created:", output_path)