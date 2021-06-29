# Wibisana2021

Imaging and single cell sequencinge analyses of super-enhancer activation mediated by NF-κB in B cells (Wibisana et al., 2021)

## Introduction

This repository contains the source codes for the sequence analysis used in the above paper.

Relative paths under "src" are shown

## scRNA-seq analysis (1_RNAseq)

Preprocessing data and alignment

- 1-1_preprocessing_alignment
  - 1_trimming_qc_alignment.sh

Downstream analysis (Clustering, pseudodose analysis, etc.)

- 1-2_singlecell_analysis
  - 1_scrnaseq_analysis.R
  - 2_pseudotime.R
  - 3_Fano_factor_calculations.R
  - mito_genes.R (Function to retrieve mitochondrial gene annotations)

## scATAC-seq analysis (2_ATACseq)

### Preprocessing of GTF file, FastQ data sampling and cellranger-atac pipeline

- 2-1_preprocessing_alignment_pipeline
  - 1_preprocessing_pipeline.sh
  - ggallus.config (cellranger-atac configuration file)
- 2-2_peak_calling_manipulation
  - 1_peak_calling.sh
  - 2_merge_bed_annotate.sh
- 2-3 SE_analysis
  - 1_SE_analysis.R

### Co-accessibility analysis

- 2-4_coaccessibility analysis
  - 0_1_chr_length.sh (get chromosome length from fasta file)
  - 1_assign_fragments.R
  - 2_cicero_main.R
  - 3_cicero_analyze_change.R
  - 4_annotate_conns.sh

### Motif analysis

- 2-5_motif_analysis
  - 1_homer_to_fimo.R (Convert Homer motif files to Fimo motif files)
  - 2_FIMO_motif_analysis.sh

## scATAC-seq and scRNA-seq combined analysis (3_RNA_atac)

- 1_GO_analysis.R
- 2_Fano_analysis.R

## Foci fitting (4_foci)

- 1_logistic_regression.R
- 2_hill_fitting.R

## Data files

### RelA-GFP foci data (../data/foci)

- 1_foci_quant.csv (All quantified foci across all dose and time points)
- 2_foci_20min.csv (Quantified foci at 20 minutes across all dose points)
- 3_median_fitting.csv (Median quantified foci at 20 minutes across all dose points for fitting)

For both files, the first column `dose` contains information of the anti-IgM concentration. Column `CELL` contains unique cell identifier per viewpoint. `N_Q570` is the target RNA spot number, while `N_Q670` is GAPDH. Column `GFP ` contains RelA-GFP foci number measured using the same method with RNA-FISH spot quantification.

### RNA-FISH data (../data/RNA_FISH)

- CD83.csv
- NFKBIA.csv

For both files, the first column `dose` contains information of the anti-IgM concentration. Column `CELL` contains unique cell identifier per viewpoint. `N_Q570` is the target RNA spot number, while `N_Q670` is GAPDH. Column `GFP ` contains RelA-GFP foci number measured using the same method with RNA-FISH spot quantification.
