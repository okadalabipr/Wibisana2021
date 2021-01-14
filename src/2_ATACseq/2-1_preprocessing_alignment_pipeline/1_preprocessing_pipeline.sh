#!/bin/sh

# Specify directories
dir_fasta="directory of genome fasta file/xxx.fa"
dir_GTF="directory of annotation GTF file/xxx.gtf"
dir_GTF_gene_name="directory of row extracted annotation GTF file/xxx_gene_name.gtf"
dir_fastq="directory containing paired-end sequence fastq files"

# Specify CPU and memory to use
N_thread="number of threads"
N_mem="GB of memory"

# Extract row with gene name only
grep "gene_name" ${dir_GTF} > Gallus_gallus.GRCg6a.96_gene_name.gtf

# Trim fastq file to 250 million reads
zcat ${dir_fastq}/xxx.fastq | head -n 1,000,000,000 > ${dir_fastq}/xxx_250M.fastq 

# Make reference file from fasta and annotation (GTF)
# config file used is available in the same directory as this script
cellranger-atac mkref NAME_OF_GENOME --config ggallus.config

# Run 10x cellranger-atac pipeline
cellranger-atac count --id=run_id --fastqs=${dir_fastq} --reference=${dir_GTF_gene_name} --localmem=${N_mem}

# by default, cellranger-atac output will be located in the home directory
