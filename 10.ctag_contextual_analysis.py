import pandas as pd


input_file = "output/final_tri_analysis.xlsx"
output_file = "output/final_CTAg_analysis.xlsx"

df = pd.read_excel(input_file)

df.rename(columns={df.columns[0]: "TRINUC"}, inplace=True)

print("Columns:", df.columns.tolist())
print("Sample data:\n", df.head())

def safe_div(a, b):
    return a / b if (b is not None and b != 0) else None

# =========================================================
#  SHEET 1: TERMINATION CODONS (TAG, TAA, TGA)
# =========================================================

stop_codons = ["TAG", "TAA", "TGA"]
df_stop = df[df["TRINUC"].isin(stop_codons)].copy()

termination_df = df_stop[[
    "TRINUC",
    "Genome_OE",
    "Coding_OE",
    "Term_OE"
]].copy()

termination_df.columns = ["Codon", "Genome", "Coding", "Termination"]

termination_sheet = termination_df.pivot_table(
    index=None,
    columns="Codon",
    values=["Genome", "Coding", "Termination"]
)

termination_sheet.columns = [
    f"{region}_{codon}" for region, codon in termination_sheet.columns
]

termination_sheet.reset_index(drop=True, inplace=True)


def get_val(motif, col):
    sub = df[df["TRINUC"] == motif]
    if sub.empty:
        return None
    return sub.iloc[0][col]

# =========================================================
#  SHEET 2: CTAG RATIOS
#  (⚠ Only works if CTAG present — else safe None)
# =========================================================

ctag_ratios = pd.DataFrame({
    "Genome_CTGA/CTAG": [safe_div(get_val("CTGA", "Genome_OE"), get_val("CTAG", "Genome_OE"))],
    "Genome_CTAA/CTAG": [safe_div(get_val("CTAA", "Genome_OE"), get_val("CTAG", "Genome_OE"))],

    "Coding_CTGA/CTAG": [safe_div(get_val("CTGA", "Coding_OE"), get_val("CTAG", "Coding_OE"))],
    "Coding_CTAA/CTAG": [safe_div(get_val("CTAA", "Coding_OE"), get_val("CTAG", "Coding_OE"))],

    "Termination_CTGA/CTAG": [safe_div(get_val("CTGA", "Term_OE"), get_val("CTAG", "Term_OE"))],
    "Termination_CTAA/CTAG": [safe_div(get_val("CTAA", "Term_OE"), get_val("CTAG", "Term_OE"))]
})

# =========================================================  
#  SHEET 3: CTAH (CTAG, CTAA, CTAC, CTAT)
# =========================================================  

ctah_codons = ["CTAG", "CTAA", "CTAC", "CTAT"]
ctah_df = df[df["TRINUC"].isin(ctah_codons)].copy()

if not ctah_df.empty:
    ctah_df = ctah_df[["TRINUC", "Genome_OE", "Coding_OE", "Term_OE"]]
    ctah_df.columns = ["Codon", "Genome", "Coding", "Termination"]

    ctah_sheet = ctah_df.pivot_table(
        index=None,
        columns="Codon",
        values=["Genome", "Coding", "Termination"]
    )

    ctah_sheet.columns = [f"{r}_{c}" for r, c in ctah_sheet.columns]
    ctah_sheet.reset_index(drop=True, inplace=True)
else:
    ctah_sheet = pd.DataFrame()

# =========================================================
#  SHEET 4: DTAG (CTAG, GTAG, ATAG, TTAG)
# =========================================================

dtag_codons = ["CTAG", "GTAG", "ATAG", "TTAG"]
dtag_df = df[df["TRINUC"].isin(dtag_codons)].copy()

if not dtag_df.empty:
    dtag_df = dtag_df[["TRINUC", "Genome_OE", "Coding_OE", "Term_OE"]]
    dtag_df.columns = ["Codon", "Genome", "Coding", "Termination"]

    dtag_sheet = dtag_df.pivot_table(
        index=None,
        columns="Codon",
        values=["Genome", "Coding", "Termination"]
    )

    dtag_sheet.columns = [f"{r}_{c}" for r, c in dtag_sheet.columns]
    dtag_sheet.reset_index(drop=True, inplace=True)
else:
    dtag_sheet = pd.DataFrame()


with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
    termination_sheet.to_excel(writer, sheet_name="Termination Codon", index=False)
    ctag_ratios.to_excel(writer, sheet_name="CTAG Ratios", index=False)

    if not ctah_sheet.empty:
        ctah_sheet.to_excel(writer, sheet_name="CTAH", index=False)

    if not dtag_sheet.empty:
        dtag_sheet.to_excel(writer, sheet_name="DTAG", index=False)

print("SUCCESS: File created at", output_file)