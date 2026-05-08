import pandas as pd
import os
import time
from Bio import Entrez

# Always provide your email when using NCBI Entrez
Entrez.email = "your_email@example.com"


data = pd.read_csv("requirements/bacterial_genome.csv")
os.makedirs("data/raw/genomes", exist_ok=True)
os.makedirs("data/raw/annotations", exist_ok=True)

print("Starting genome downloads...\n")

for index, row in data.iterrows():

    accession = row["accession"]
    organism = row["organism"]

    print(f"Downloading: {organism} ({accession})")

    try:
        # -------- Download FASTA genome -------- 
        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="fasta",
            retmode="text"
        )

        fasta_data = handle.read()

        with open(f"data/raw/genomes/{organism}.fna", "w") as f:
            f.write(fasta_data)

        handle.close()

        # -------- Download GenBank annotation --------
        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="gbwithparts",
            retmode="text"
        )

        gb_data = handle.read()

        with open(f"data/raw/annotations/{organism}.gb", "w") as f:
            f.write(gb_data)

        handle.close()

        print("Download completed\n")

        time.sleep(1)
    except Exception as e:
        print("Error downloading", organism, ":", e)

print("All genomes downloaded successfully.")