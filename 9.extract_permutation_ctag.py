import pandas as pd
import os

input_file = "output/final_CTAg_analysis.xlsx"
output_file = "output/CTAG_permutation_final.xlsx"

# -------- CTAG PERMUTATIONS --------
ctag_permutations = [
"CTAG","ACGT","ACTG","AGCT","AGTC",
"ATCG","ATGC","CAGT","CATG","CGAT",
"CGTA","CTGA","GACT","GATC","GCAT",
"GCTA","GTAC","GTCA","TACG","TAGC",
"TCAG","TCGA","TGAC","TGCA"
]


xls = pd.ExcelFile(input_file)

genome_rows = []
coding_rows = []
termination_rows = []

for sheet in xls.sheet_names:

    df = pd.read_excel(xls, sheet_name=sheet)

    organism = sheet

    df = df[df["Tetranucleotide"].isin(ctag_permutations)]

    # ---- GENOME ----
    genome_dict = {"BACTERIAL SPECIES": organism}
    for _, row in df.iterrows():
        genome_dict[row["Tetranucleotide"]] = row["Genome_OE"]
    genome_rows.append(genome_dict)

    # ---- CODING ----
    coding_dict = {"BACTERIAL SPECIES": organism}
    for _, row in df.iterrows():
        coding_dict[row["Tetranucleotide"]] = row["Coding_OE"]
    coding_rows.append(coding_dict)

    # ---- TERMINATION ----
    term_dict = {"BACTERIAL SPECIES": organism}
    for _, row in df.iterrows():
        term_dict[row["Tetranucleotide"]] = row["Term_OE"]
    termination_rows.append(term_dict)


genome_df = pd.DataFrame(genome_rows)
coding_df = pd.DataFrame(coding_rows)
termination_df = pd.DataFrame(termination_rows)


columns_order = ["BACTERIAL SPECIES"] + ctag_permutations

genome_df = genome_df[columns_order]
coding_df = coding_df[columns_order]
termination_df = termination_df[columns_order]

with pd.ExcelWriter(output_file, engine="openpyxl") as writer:

    genome_df.to_excel(writer, sheet_name="genome_permutation", index=False)
    coding_df.to_excel(writer, sheet_name="coding_permutation", index=False)
    termination_df.to_excel(writer, sheet_name="termination_permutation", index=False)

print("\n DONE! Final permutation file created:", output_file)