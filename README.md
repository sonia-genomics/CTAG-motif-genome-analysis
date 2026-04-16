# CTAG Motif Genome Analysis

Analysis of **CTAG under-representation in bacterial genomes** using a multi-level bioinformatics pipeline integrating genome-wide, coding, and termination site analyses.

---

## Overview

DNA sequences in bacterial genomes exhibit **non-random compositional patterns**. One such phenomenon is the consistent **under-representation of the CTAG tetranucleotide**.

This project presents a computational pipeline to investigate the **distribution, frequency, and biological significance of CTAG motifs** across bacterial genomes.

---

## Objectives

* Quantify **Observed/Expected (O/E) ratios** of tetranucleotides
* Investigate **CTAG under-representation** across genomes
* Compare motif distribution across:

  * Whole genome
  * Coding regions
  * Termination sites
* Analyze **contextual and positional effects** influencing CTAG frequency

---

## Methodology

### 1. Genome Processing

* Genome sequences loaded and standardized
* Quality control performed on input datasets

### 2. Motif Frequency Analysis

* Generated all possible **256 tetranucleotides**

* Computed:

  * Observed counts (sliding window approach)
  * Expected frequencies based on base composition

* O/E ratio:

```
O/E = Observed Frequency / Expected Frequency
```

---

### 3. Region-Specific Analysis

Motif frequencies were computed across:

* **Genome-wide sequences**
* **Coding regions** (from GenBank annotations)
* **Termination sites** (derived from experimental datasets)

---

### 4. CTAG-Specific Analysis

* Contextual analysis of CTAG and related motifs

* Comparative analysis with:

  * **CTAH group:** CTAG, CTAA, CTAC, CTAT
  * **DTAG group:** CTAG, GTAG, ATAG, TTAG

* Ratio-based metrics:

  * (CTAA + CTGA) / CTAG
  * Compared with stop codon ratios (TAA + TGA) / TAG

---

### 5. Statistical & Comparative Analysis

* Comparative analysis across **23 bacterial genomes**
* Identification of:

  * Under-represented motifs (O/E < 1)
  * Region-specific enrichment patterns

---

### 6. Visualization

* Generated **boxplots** comparing:

  * Genome vs Coding vs Termination regions
* Visualization of motif distribution and ratio differences

---

## Key Findings

* **CTAG is consistently under-represented** across bacterial genomes (O/E < 1)
* CTAG frequently shows the **lowest O/E value among tetranucleotides**
* The submotif **TAG (stop codon)** is also under-represented
* Addition of **5′ cytosine (forming CTAG)** further reduces frequency
* CTAG distribution is **highly heterogeneous**, with localized enrichment regions
* Base composition plays a key role in **motif inhomogeneity**

---

## Biological Significance

These findings suggest that CTAG under-representation is **non-random and biologically driven**, potentially associated with:

* **Transcription termination mechanisms**
* **Genomic stability constraints**
* **Regulatory sequence selection in bacterial genomes**

---

## Technologies Used

* Python
* Biopython
* Pandas, NumPy
* Matplotlib
* Linux

---

## Pipeline Modules

| Step                 | Script                                  | Function             |
| -------------------- | --------------------------------------- | -------------------- |
| Data Acquisition     | `download_data.py`                      | Genome retrieval     |
| Quality Control      | `genome_qc.py`                          | Data filtering       |
| Gene Processing      | `extract_genes_from_genbank.py`         | Gene extraction      |
| Coordinates          | `gene_coordinates.py`                   | Gene mapping         |
| Motif Analysis       | `calculate_obs_exp_tetra.py`            | O/E computation      |
| Context Analysis     | `ctag_contextual_analysis.py`           | CTAG-specific ratios |
| Termination Analysis | `ctag_termination_analysis_pipeline.py` | Comparative analysis |
| Visualization        | `visualization.py`                      | Plot generation      |

---

## Future Directions

* Extend analysis to **eukaryotic genomes**
* Integrate **multi-omics datasets**
* Apply findings to **disease-associated genomic regions**
* Develop scalable workflows using **Snakemake or Nextflow**

---

## Author

**Sonia**
MSc Biotechnology | Bioinformatics & Genomics
Aspiring PhD researcher in **genomics and genetic diseases**
