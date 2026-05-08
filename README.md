# CTAG Motif Genome Analysis

## Analysis of CTAG Under-Representation, Local Enrichment, and Genome Organization in Bacterial Genomes

A large-scale computational genomics pipeline for investigating the biological significance, suppression patterns, positional distribution, and local clustering behavior of the CTAG tetranucleotide motif across bacterial genomes.

This project integrates comparative genomics, motif frequency analysis, statistical genomics, genome architecture analysis, entropy measurements, phylogenetic comparisons, and reproducible workflow automation using Python and Snakemake.

---

# Overview

DNA sequences in bacterial genomes exhibit strong non-random compositional patterns. One such conserved phenomenon is the significant under-representation of the CTAG tetranucleotide motif across diverse bacterial species.

This project presents a multi-level computational framework to investigate:

- CTAG suppression across bacterial genomes
- Contextual sequence effects surrounding CTAG
- Coding vs non-coding motif distribution
- Termination codon-associated motif bias
- Local CTAG clustering and hotspot regions
- Base heterogeneity associated with CTAG-rich regions
- Genome-wide motif architecture and entropy
- Comparative genome organization patterns

The project expands beyond simple motif counting and explores the relationship between CTAG distribution, genome organization, sequence heterogeneity, and possible biological constraints acting on bacterial genomes.

---

# Highlights

- Multi-genome bacterial comparative genomics pipeline
- Automated reproducible workflow using Snakemake
- Genome-wide tetranucleotide and trinucleotide analysis
- Observed/Expected (O/E) motif frequency calculations
- CTAG hotspot and local clustering detection
- Base heterogeneity and entropy analysis
- CTAG vs GC-content correlation analysis
- PCA and hierarchical genome clustering
- Phylogenetic comparative analysis
- Publication-style scientific visualizations
- Functional annotation of CTAG-enriched regions
- Comparative motif architecture across coding and termination regions

---

# Biological Motivation

The CTAG tetranucleotide is one of the most strongly under-represented motifs in many bacterial genomes.

Interestingly:

- TAG itself functions as a termination codon
- Addition of a 5′ cytosine (forming CTAG) further suppresses motif frequency
- CTAG motifs often appear sparsely distributed genome-wide
- However, specific genomic regions show unexpectedly high local CTAG abundance

These localized CTAG-enriched regions frequently overlap with regions exhibiting:

- High base heterogeneity
- Irregular nucleotide composition
- Distinct genomic architecture

Such findings suggest that CTAG suppression and clustering may be biologically driven rather than random genomic events.

---

# Objectives

## Primary Objectives

- Quantify CTAG under-representation across bacterial genomes
- Compute tetranucleotide and trinucleotide O/E ratios
- Compare motif distributions across:
  - Whole genomes
  - Coding sequences
  - Termination regions
- Investigate positional and contextual effects associated with CTAG

---

## Advanced Genomic Objectives

- Detect local CTAG clustering regions
- Identify CTAG hotspot enrichment patterns
- Analyze genome-wide sequence heterogeneity
- Measure genomic entropy and compositional variability
- Study correlation between CTAG density and GC content
- Perform phylogenetic comparison of CTAG organization
- Cluster bacterial genomes based on motif architecture

---

# Project Structure

```text
CTAG-motif-genome-analysis/
│
├── scripts/
│
├── workflow/
│   └── Snakefile
│
├── requirements/
│   └── bacterial_genome.csv
│
├── output/
│
├── results/
│   ├── figures/
│   ├── tables/
│   └── reports/
│
├── docs/
│
├── environment.yml
├── README.md
└── run_pipeline.sh
```

---

# Pipeline Workflow

```text
Genome Download
        ↓
Genome Quality Control
        ↓
Gene Coordinate Extraction
        ↓
Gene Sequence Extraction
        ↓
Coding Sequence Construction
        ↓
Termination Site Extraction
        ↓
Tetranucleotide Analysis
        ↓
Trinucleotide Analysis
        ↓
CTAG Contextual Analysis
        ↓
Statistical Analysis
        ↓
Heatmaps + PCA + Clustering
        ↓
Local CTAG Clustering
        ↓
Base Heterogeneity Analysis
        ↓
Entropy & GC Correlation Analysis
        ↓
Phylogenetic Comparative Analysis
        ↓
Hotspot Annotation
        ↓
Comparative Genome Summary
        ↓
Final Automated Report
```

---

# Pipeline Modules

| Step | Script | Purpose |
|---|---|---|
| 01 | `01_download_data.py` | Genome download from NCBI |
| 02 | `02_genome_qc.py` | Genome quality control |
| 03 | `03_gene_coordinates.py` | CDS coordinate extraction |
| 04 | `04_extract_genes_from_genbank.py` | Gene sequence extraction |
| 05 | `05_combine_gene_seq.py` | Combined coding FASTA generation |
| 06 | `06_extract_termination_sites.py` | Termination motif extraction |
| 07 | `07_concat_fasta_seq.py` | Combined FASTA generation |
| 08 | `08_calculate_obs_exp_tetra.py` | Tetranucleotide O/E analysis |
| 09 | `09_calculate_obs_exp_tri.py` | Trinucleotide O/E analysis |
| 10 | `10_extract_permutation_ctag.py` | CTAG permutation analysis |
| 11 | `11_ctag_termination_analysis_pipeline.py` | Comparative codon analysis |
| 12 | `12_visualization.py` | Publication-style figures |
| 13 | `13_statistical_analysis.py` | Statistical genomics analysis |
| 14 | `14_heatmap_analysis.py` | Genome heatmap generation |
| 15 | `15_pca_analysis.py` | PCA clustering analysis |
| 16 | `16_genome_clustering.py` | Hierarchical clustering |
| 17 | `17_ctag_publication_visualization.py` | Advanced scientific plots |
| 18 | `18_ctag_local_clustering_analysis.py` | Local CTAG enrichment analysis |
| 19 | `19_base_heterogeneity_analysis.py` | Base heterogeneity analysis |
| 20 | `20_ctag_vs_gc_correlation.py` | GC vs CTAG correlation |
| 21 | `21_genome_entropy_analysis.py` | Genome entropy analysis |
| 22 | `22_phylogenetic_ctag_analysis.py` | Phylogenetic CTAG comparison |
| 23 | `23_ctag_hotspot_annotation.py` | CTAG hotspot annotation |
| 24 | `24_ctag_genome_browser_tracks.py` | Genome browser track generation |
| 25 | `25_ctag_hotspot_functional_enrichment.py` | Functional hotspot analysis |
| 26 | `26_ctag_comparative_summary.py` | Comparative genomics summary |
| 27 | `27_generate_final_report.py` | Automated report generation |

---

# Statistical Analyses Performed

- Observed/Expected motif analysis
- Genome-wide motif distribution analysis
- CTAG enrichment detection
- Sliding window motif analysis
- PCA dimensionality reduction
- Hierarchical clustering
- Entropy measurements
- Base heterogeneity calculations
- Correlation analysis
- Comparative genomics statistics

---

# Key Findings

- CTAG is consistently under-represented across bacterial genomes
- CTAG often exhibits the lowest O/E ratio among tetranucleotides
- TAG trinucleotide suppression contributes to CTAG depletion
- Addition of a 5′ cytosine significantly enhances motif suppression
- CTAG motifs are unevenly distributed across genomes
- Specific genomic regions show strong local CTAG enrichment
- CTAG hotspot regions frequently overlap with high base heterogeneity
- Genome organization influences motif architecture and suppression patterns

---

# Scientific Significance

This project demonstrates that CTAG suppression is likely a biologically constrained genomic phenomenon rather than a random sequence property.

The observed relationship between:

- CTAG clustering
- Base heterogeneity
- Genomic entropy
- Sequence organization

suggests that bacterial genomes may maintain selective pressures influencing motif architecture and local nucleotide composition.

---

# Reproducibility

The entire workflow is fully automated using Snakemake.

All analyses are reproducible from a single genome accession input file:

```text
requirements/bacterial_genome.csv
```

Run the complete pipeline:

```bash
snakemake --cores 4
```

---

# Installation

## Clone Repository

```bash
git clone https://github.com/sonia-genomics/CTAG-motif-genome-analysis.git

cd CTAG-motif-genome-analysis
```

---

## Create Conda Environment

```bash
conda env create -f environment.yml

conda activate ctag-analysis
```

---

# Required Python Packages

- Python
- Biopython
- Pandas
- NumPy
- SciPy
- Matplotlib
- Seaborn
- Scikit-learn
- OpenPyXL
- Snakemake

---

# Input File

```text
requirements/bacterial_genome.csv
```

Example:

| id | Organism | family | accession |genome_size |

|---|---|---|---|---|
| 1 | Escherichia coli str. K-12 | Enterobacteriaceae | NC_000913.3 | 4641652 |
| 2 | Anabaena cylindrica PCC 7122 | Nostocaceae | NC_019771.1 | 6395836 |


---

# Example Outputs

## Tables

- Genome QC statistics
- Tetranucleotide O/E tables
- CTAG comparative statistics
- Entropy analysis tables
- Phylogenetic summary tables

---

## Figures

- CTAG boxplots
- Heatmaps
- PCA plots
- Genome dendrograms
- CTAG hotspot maps
- Base heterogeneity plots
- GC correlation figures

---

# Future Directions

- Expansion to archaeal and eukaryotic genomes
- Comparative analysis across pathogenic bacteria
- Integration with transcriptomics datasets
- Long-read genome architecture analysis
- Evolutionary modeling of motif suppression
- Genome-scale motif interaction networks
- Large-scale comparative microbial genomics

---

# Computational Skills Demonstrated

- Comparative genomics
- Statistical genomics
- Sequence analysis
- Workflow automation
- Reproducible research
- Genome-scale data processing
- Bioinformatics pipeline development
- Scientific visualization
- Computational biology research design

---

# Author

## Sonia

MSc Biotechnology  

Open to:

Research collaborations
Bioinformatics and genomics roles
Computational biology opportunities
Fully funded PhD positions