# Publication-Style CTAG Cluster Visualization
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import os

input_file = "output/CTAG_permutation_final.xlsx"
output_dir = "results/figures"

os.makedirs(output_dir, exist_ok=True)

# Use genome permutation sheet
df = pd.read_excel(
    input_file,
    sheet_name="genome_permutation"
)


# CTAG PERMUTATIONS
motifs = [
    "CTAG","ACGT","ACTG","AGCT","AGTC",
    "ATCG","ATGC","CAGT","CATG","CGAT",
    "CGTA","CTGA","GACT","GATC","GCAT",
    "GCTA","GTAC","GTCA","TACG","TAGC",
    "TCAG","TCGA","TGAC","TGCA"
]


# FIGURE LAYOUT
n_species = len(df)

ncols = 4
nrows = math.ceil(n_species / ncols)

fig, axes = plt.subplots(
    nrows,
    ncols,
    figsize=(20, 4 * nrows)
)

axes = axes.flatten()

# COLORS

observed_color = "royalblue"
expected_color = "darkorange"

# PLOT EACH ORGANISM

for idx, row in df.iterrows():

    ax = axes[idx]

    organism = row["BACTERIAL SPECIES"]

    values = []

    for motif in motifs:

        if motif in row:
            values.append(row[motif])
        else:
            values.append(np.nan)

    # EXPECTED LINE
    expected = sorted(values, reverse=True)

    # OBSERVED LINE
    observed = values

    # PLOT
    ax.plot(
        range(len(motifs)),
        observed,
        marker='o',
        linewidth=1.5,
        label='Observed CTAG',
        color=observed_color
    )

    ax.plot(
        range(len(motifs)),
        expected,
        marker='o',
        linewidth=1.5,
        label='Expected CTAG',
        color=expected_color
    )

    # HIGHLIGHT CTAG
    try:
        ctag_index = motifs.index("CTAG")

        ax.scatter(
            ctag_index,
            observed[ctag_index],
            s=80,
            edgecolor='black',
            linewidth=1,
            zorder=5
        )

        ax.annotate(
            'CTAG',
            (ctag_index, observed[ctag_index]),
            textcoords="offset points",
            xytext=(0, 10),
            ha='center',
            fontsize=7
        )

    except:
        pass

    ax.set_title(
        organism,
        fontsize=8,
        fontweight='bold'
    )

    ax.set_xticks(range(len(motifs)))
    ax.set_xticklabels(
        motifs,
        rotation=90,
        fontsize=5
    )

    ax.set_ylabel("O/E", fontsize=7)

    ax.grid(True, linestyle='--', alpha=0.4)

# REMOVE EMPTY AXES
for j in range(idx + 1, len(axes)):
    fig.delaxes(axes[j])

# GLOBAL TITLE
fig.suptitle(
    "CTAG Local Over-Representation in Bacterial Genomes",
    fontsize=18,
    fontweight='bold'
)

# GLOBAL LEGEND
handles, labels = axes[0].get_legend_handles_labels()

fig.legend(
    handles,
    labels,
    loc='lower center',
    ncol=2,
    fontsize=10
)

# LAYOUT
plt.tight_layout(rect=[0, 0.05, 1, 0.96])

output_file = os.path.join(
    output_dir,
    "ctag_cluster_publication_style.png"
)

plt.savefig(
    output_file,
    dpi=600,
    bbox_inches='tight'
)

plt.close()

print("SUCCESS: Figure saved")
print(output_file)
