#!/bin/bash
#=================================================
echo "01_download_data.py"
#=================================================
python scripts/01_download_data.py
#=================================================
echo "02_genome_qc.py"
#=================================================
python scripts/02_genome_qc.py
#=================================================
echo "03_gene_coordinates.py"
#=================================================
python scripts/03_gene_coordinates.py
#=================================================
echo "04_extract_genes_from_genbank.py"
#=================================================
python scripts/04_extract_genes_from_genbank.py
#=================================================
echo "05_combine_gene_seq.py"
#=================================================
python scripts/05_combine_gene_seq.py
#=================================================
echo "06_extract_termination_sites.py"
#=================================================
python scripts/06_extract_termination_sites.py
#=================================================
echo "07_concat_fasta_seq.py"
#=================================================
python scripts/07_concat_fasta_seq.py
#=================================================
echo "08_calculate_obs_exp_tetra.py"
#=================================================
python scripts/08_calculate_obs_exp_tetra.py
#=================================================
echo "09_calculate_obs_exp_tri.py"
#=================================================
python scripts/09_calculate_obs_exp_tri.py
#=================================================
echo "10_extract_permutation_ctag.py"
#=================================================
python scripts/10_extract_permutation_ctag.py
#=================================================
echo "11_ctag_termination_analysis_pipeline.py"
#=================================================
python scripts/11_ctag_termination_analysis_pipeline.py
#=================================================
echo "12_visualization.py"
#=================================================
python scripts/12_visualization.py
#=================================================
echo "13_statistical_analysis.py"
#=================================================
python scripts/13_statistical_analysis.py
#=================================================
echo "14_heatmap_analysis.py"
#=================================================
python scripts/14_heatmap_analysis.py
#=================================================
echo "15_pca_analysis.py"
#=================================================
python scripts/15_pca_analysis.py
#=================================================
echo "16_genome_clustering.py"
#=================================================
python scripts/16_genome_clustering.py
#=================================================
echo "17_ctag_publication_visualization.py"
#=================================================
python scripts/17_ctag_publication_visualization.py
#=================================================
echo "18_ctag_local_clustering_analysis.py"
#=================================================
python scripts/18_ctag_local_clustering_analysis.py
#=================================================
echo "19_base_heterogeneity_analysis.py"
#=================================================
python scripts/19_base_heterogeneity_analysis.py
#=================================================
echo "20_ctag_vs_gc_correlation.py"
#=================================================
python scripts/20_ctag_vs_gc_correlation.py
#=================================================
echo "21_genome_entropy_analysis.py"
#=================================================
python scripts/21_genome_entropy_analysis.py
#=================================================
echo "22_phylogenetic_ctag_analysis.py"
#=================================================
python scripts/22_phylogenetic_ctag_analysis.py
#=================================================
echo "23_ctag_hotspot_annotation.py"
#=================================================
python scripts/23_ctag_hotspot_annotation.py
#=================================================
echo "24_ctag_genome_browser_tracks.py"
#=================================================
python scripts/24_ctag_genome_browser_tracks.py
#=================================================
echo "25_ctag_hotspot_functional_enrichment.py"
#=================================================
python scripts/25_ctag_hotspot_functional_enrichment.py
#=================================================
echo "26_ctag_comparative_summary.py"
#=================================================
python scripts/26_ctag_comparative_summary.py
#=================================================
echo "27_generate_final_report.py"
#=================================================
python scripts/27_generate_final_report.py