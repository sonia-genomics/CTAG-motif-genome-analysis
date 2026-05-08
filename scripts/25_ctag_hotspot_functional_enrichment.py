
# Extract functional annotations associated
# with CTAG hotspot regions.
# =========================================================

import os
import pandas as pd
from collections import Counter

hotspot_dir = "output/ctag_hotspots"

output_dir = "output/functional_enrichment"
os.makedirs(output_dir, exist_ok=True)

summary_output = os.path.join(
    output_dir,
    "functional_category_summary.csv"
)

# FUNCTIONAL KEYWORDS
categories = {
    "Transcription": [
        "transcription",
        "sigma",
        "rna polymerase"
    ],

    "Replication": [
        "dna polymerase",
        "replication",
        "helicase"
    ],

    "Membrane": [
        "membrane",
        "transport",
        "channel"
    ],

    "Metabolism": [
        "metabolism",
        "synthase",
        "dehydrogenase"
    ],

    "Mobile_Elements": [
        "transposase",
        "integrase",
        "phage"
    ]
}

# ANALYZE HOTSPOT FILES
all_results = []

for file in os.listdir(hotspot_dir):

    if not file.endswith(".csv"):
        continue

    organism = file.replace("_ctag_hotspots.csv", "")

    path = os.path.join(hotspot_dir, file)

    df = pd.read_csv(path)

    functions = []

    for genes in df["Nearby_Genes"].dropna():

        text = str(genes).lower()

        matched = False

        for category, keywords in categories.items():

            for keyword in keywords:

                if keyword in text:

                    functions.append(category)
                    matched = True
                    break

        if not matched:
            functions.append("Other")

    counts = Counter(functions)

    for category, count in counts.items():

        all_results.append({
            "Organism": organism,
            "Category": category,
            "Count": count
        })


result_df = pd.DataFrame(all_results)

result_df.to_csv(summary_output, index=False)

print("\n=================================================")
print("FUNCTIONAL ENRICHMENT ANALYSIS COMPLETED")
print("=================================================")