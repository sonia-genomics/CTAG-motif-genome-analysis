import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

input_file = "output/final_analysis.xlsx"
output_dir = os.path.join(os.path.dirname(input_file), "thesis_figures")
os.makedirs(output_dir, exist_ok=True)


t4 = pd.read_excel(input_file, sheet_name="Table 4")
t5 = pd.read_excel(input_file, sheet_name="Table 5")
t6 = pd.read_excel(input_file, sheet_name="Table 6")
t7 = pd.read_excel(input_file, sheet_name="Table 7")


def clean_data(df):
    data = df.iloc[:, 1:]
    return data.apply(pd.to_numeric, errors='coerce')

def prepare(df, codons, region):
    df.columns = codons
    melted = df.melt(var_name="Codon", value_name="OE")
    melted["Region"] = region
    return melted.dropna()

# =========================================================
# COLOR MAP (CODON-BASED)
# =========================================================
COLOR_MAP = {
    "TAG": "blue",
    "CTAG": "blue",

    "TAA": "red",
    "CTAA": "red",

    "TGA": "green",
    "CTGA": "green",

    "GTAG": "yellow",
    "ATAG": "green",
    "TTAG": "red"
}


def plot_seaborn(data, codons, title, filename):

    n = len(codons)

    genome = data.iloc[:, 0:n]
    coding = data.iloc[:, n:2*n]
    term = data.iloc[:, 2*n:3*n]

    df_all = pd.concat([
        prepare(genome, codons, "GENOME"),
        prepare(coding, codons, "CODING"),
        prepare(term, codons, "TERMINATION")
    ])

    g = sns.catplot(
        data=df_all,
        x="Codon",
        y="OE",
        hue="Codon",
        col="Region",
        kind="box",
        palette=COLOR_MAP,
        height=5,
        aspect=1,
        sharey=False,
        showfliers=False
        
    )

    g.fig.suptitle(title, y=1.05)

    for ax in g.axes.flat:
        ax.axhline(1, linestyle="--", color="orange")
        ax.set_ylabel("O/E VALUE")
        ax.set_xlabel("Codon")

    plt.tight_layout()

    save_path = os.path.join(output_dir, filename)
    plt.savefig(save_path, dpi=300)
    plt.close()

    print("Seaborn saved:", save_path)


# =========================================================
# GENERATE FIGURES
# =========================================================

# Fig 3
data = clean_data(t4)
plot_seaborn(data, ["TAG","TAA","TGA"], "Fig 3: Termination Codons", "Fig3_seaborn.png")

# Fig 4
data = clean_data(t5)
plot_seaborn(data, ["CTAG","CTAA","CTGA"], "Fig 4: Cytosine Preceding", "Fig4_seaborn.png")

# Fig 5
data = clean_data(t7)
plot_seaborn(data, ["CTAG","GTAG","ATAG","TTAG"], "Fig 5: CTAG vs DTAG", "Fig5_seaborn.png")

print("\n ALL FIGURES GENERATED SUCCESSFULLY!")