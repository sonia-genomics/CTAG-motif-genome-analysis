import os
import pandas as pd
from Bio import SeqIO
from collections import Counter
import itertools

# -------- PATHS --------
genome_folder = "output/genomes_singleline"
genes_folder = "output/genes"
termination_file = "output/last4_nt.xlsx"

output_path = "output/final_CTAg_analysis.xlsx"

bases = ["A","T","G","C"]
tetranucleotides = ["".join(p) for p in itertools.product(bases, repeat=4)]


termination_xls = pd.ExcelFile(termination_file)


with pd.ExcelWriter(output_path, engine="openpyxl") as writer:

    for genome_file in os.listdir(genome_folder):

        if genome_file.endswith(".fna"):

            organism = genome_file.replace(".fna","")
            genome_path = os.path.join(genome_folder, genome_file)

            print("Processing:", organism)

            # =========================
            # GENOME
            # =========================
            record = SeqIO.read(genome_path, "fasta")
            genome_seq = str(record.seq).upper()
            genome_len = len(genome_seq)

            kmer_counts_genome = Counter(genome_seq[i:i+4] for i in range(len(genome_seq)-3))

            # base freq
            fa = genome_seq.count("A")/genome_len
            ft = genome_seq.count("T")/genome_len
            fg = genome_seq.count("G")/genome_len
            fc = genome_seq.count("C")/genome_len

            genome_expected_dict = {}
            for tet in tetranucleotides:
                f1 = {"A":fa,"T":ft,"G":fg,"C":fc}[tet[0]]
                f2 = {"A":fa,"T":ft,"G":fg,"C":fc}[tet[1]]
                f3 = {"A":fa,"T":ft,"G":fg,"C":fc}[tet[2]]
                f4 = {"A":fa,"T":ft,"G":fg,"C":fc}[tet[3]]
                genome_expected_dict[tet] = f1*f2*f3*f4*genome_len

            # =========================
            # CODING
            # =========================
            genes_file = os.path.join(genes_folder, organism + "_genes.fasta")

            coding_seq = ""
            if os.path.exists(genes_file):
                for rec in SeqIO.parse(genes_file, "fasta"):
                    coding_seq += str(rec.seq).upper()

            coding_len = len(coding_seq)
            kmer_counts_coding = Counter(coding_seq[i:i+4] for i in range(len(coding_seq)-3))

            # =========================
            # TERMINATION (FIXED)
            # =========================
            term_counts = Counter()
            term_len = 0

            matched_sheet = None
            for sheet in termination_xls.sheet_names:
                if sheet in organism or organism in sheet:
                    matched_sheet = sheet
                    break

            if matched_sheet:
                term_df = pd.read_excel(termination_xls, sheet_name=matched_sheet)

                last4_list = term_df["Last_4_nt"].dropna().str.upper().tolist()

                term_counts = Counter(last4_list)

                term_len = len(last4_list) * 4  # total nt

            # =========================
            # FINAL TABLE
            # =========================
            data = []

            for tet in tetranucleotides:

                # GENOME
                obs_g = kmer_counts_genome.get(tet, 0)
                exp_g = genome_expected_dict.get(tet, 0)
                oe_g = obs_g / exp_g if exp_g > 0 else 0

                # CODING
                obs_c = kmer_counts_coding.get(tet, 0)
                exp_c = exp_g * (coding_len / genome_len)
                oe_c = obs_c / exp_c if exp_c > 0 else 0

                # TERMINATION
                obs_t = term_counts.get(tet, 0)
                exp_t = exp_g * (term_len / genome_len)
                oe_t = obs_t / exp_t if exp_t > 0 else 0

                data.append([
                    tet,genome_len,coding_len,term_len,
                    obs_g, exp_g, oe_g,
                    obs_c, exp_c, oe_c,
                    obs_t, exp_t, oe_t
                ])

            df = pd.DataFrame(data, columns=[
                "Tetranucleotide","genome","coding","termination",
                "Genome_obs","Genome_exp","Genome_OE",
                "Coding_obs","Coding_exp","Coding_OE",
                "Term_obs","Term_exp","Term_OE"
            ])

            sheet_name = organism[:31]
            df.to_excel(writer, sheet_name=sheet_name, index=False)

            print("Saved sheet:", sheet_name)

print("\n DONE! File saved at:", output_path)