import os
from Bio import SeqIO

genome_folder = "output/genomes"
annotation_folder = "output/annotations"
output_folder = "output/genes"

os.makedirs(output_folder, exist_ok=True)

for genome_file in os.listdir(genome_folder):

    if genome_file.endswith(".fna"):

        organism = genome_file.replace(".fna","")

        genome_path = os.path.join(genome_folder, genome_file)
        gb_path = os.path.join(annotation_folder, organism + ".gb")
        output_file = os.path.join(output_folder, organism + "_genes.fasta")

        if not os.path.exists(gb_path):
            print("GenBank file missing for:", organism)
            continue

        print("Processing:", organism)

        # -----------------------------
        # STEP 1: Load genome
        # -----------------------------

        genome_record = SeqIO.read(genome_path, "fasta")
        genome_seq = genome_record.seq

        print("Genome length:", len(genome_seq))

        # -----------------------------
        # STEP 2: Read annotation
        # -----------------------------

        gb_record = SeqIO.read(gb_path, "genbank")

        genes = []

        for feature in gb_record.features:

            if feature.type == "CDS":

                start = int(feature.location.start) + 1
                end = int(feature.location.end)
                strand = feature.location.strand

                gene_seq = genome_seq[start-1:end]

                if strand == -1:
                    gene_seq = gene_seq.reverse_complement()

                # get gene name
                if "gene" in feature.qualifiers:
                    gene_name = feature.qualifiers["gene"][0]

                elif "locus_tag" in feature.qualifiers:
                    gene_name = feature.qualifiers["locus_tag"][0]

                else:
                    gene_name = "unknown_gene"

                genes.append((gene_name, start, end, strand, gene_seq))

        print("Total CDS extracted:", len(genes))

        # -----------------------------
        # STEP 3: Save gene sequences
        # -----------------------------

        with open(output_file, "w") as out:

            for gene_name, start, end, strand, seq in genes:

                strand_symbol = "+" if strand == 1 else "-"

                out.write(f">{gene_name} {start}:{end} ({strand_symbol})\n")
                out.write(str(seq) + "\n")

        print("Gene sequences saved to:", output_file)