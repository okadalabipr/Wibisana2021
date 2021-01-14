# Wibisana2021

Wibisana, JN et al., DT40 paper

## Introduction

This repository contains the source codes for the sequence analysis used in the above paper.

Relative paths under "src" are shown

### scRNA-seq analysis (1_RNAseq)

Preprocessing data and alignment

- 1-1_preprocessing_alignment
  - 1_trimming_qc_alignment.sh

Downstream analysis (Clustering, pseudodose analysis, etc.)

- 1-2_singlecell_analysis

### scATAC-seq analysis (2_ATACseq)

Preprocessing of GTF file, FastQ data sampling and cellranger-atac pipeline

- 2-1_preprocessing_alignment_pipeline
  - 1_preprocessing_pipeline.sh
  - ggallus.config (cellranger-atac configuration file)


### 
