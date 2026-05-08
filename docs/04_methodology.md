# CHAPTER 4 — METHODOLOGY

The project implements a complete computational genomics workflow
for CTAG motif analysis.

---

# Pipeline Overview

The workflow includes:

1. Genome download
2. Quality control
3. Gene extraction
4. Termination site extraction
5. Tetranucleotide analysis
6. Trinucleotide analysis
7. Statistical analysis
8. Heatmap analysis
9. PCA analysis
10. Genome clustering
11. Local CTAG clustering
12. Base heterogeneity analysis
13. GC correlation analysis
14. Entropy analysis
15. Hotspot annotation
16. Genome browser track generation
17. Functional enrichment analysis
18. Comparative genomics summary
19. Automated report generation

---

# Observed/Expected Analysis

Tetranucleotide abundance was measured using:

O/E = Observed Frequency / Expected Frequency

where:

- Observed frequency = actual motif count
- Expected frequency = probability-based expected count

---

# Statistical Analysis

The pipeline includes:

- Z-score analysis
- Permutation analysis
- Comparative ranking
- Correlation analysis
- Clustering
- Principal Component Analysis

---

# Visualization

Generated visualizations include:

- Boxplots
- Heatmaps
- PCA plots
- Genome dendrograms
- Local clustering maps
- Entropy plots
- Browser tracks

---

# Workflow Automation

The entire pipeline is automated using Snakemake
for reproducibility and scalability.