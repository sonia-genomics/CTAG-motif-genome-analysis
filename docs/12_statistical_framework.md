# STATISTICAL FRAMEWORK

The project uses multiple statistical approaches to evaluate
CTAG suppression patterns.

---

# Observed/Expected Ratio

Primary metric:

O/E = Observed Frequency / Expected Frequency

Interpretation:

- O/E < 1 → under-representation
- O/E > 1 → over-representation

---

# Permutation Analysis

Permutation-based analysis evaluates whether CTAG suppression
is stronger than expected relative to:

- Base composition
- Random motif distribution
- Other tetranucleotide permutations

---

# Comparative Ranking

Each tetranucleotide is ranked according to O/E value.

This identifies:

- Most suppressed motifs
- Most enriched motifs
- Genome-specific motif behavior

---

# Correlation Analysis

The workflow evaluates relationships between:

- CTAG frequency
- GC content
- Genome entropy
- Base heterogeneity

---

# Clustering Analysis

Implemented analyses include:

- PCA
- Hierarchical clustering
- Distance-based genome grouping

These methods reveal similarity between bacterial genomes.

---

# Local Clustering Statistics

Poisson-based local enrichment analysis identifies:

- CTAG hotspots
- Significant motif accumulation regions
- Non-random local distributions