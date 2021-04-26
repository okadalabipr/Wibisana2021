#!/bin/sh

# specify input and output directories
dir_bed="directory containing input bed file/xxxx.bed"
dir_out="output directory"
dir_motifs="directory containing meme format motifs"

# Specify fasta (fa) directory 
dir_fasta="directory containing fasta file/xxxx.fa"

# Make fasta from bed file of regions to analyze and visualize
bedtools getfasta \
-fi ${dir_fasta} \
-bed ${dir_bed} \
-fo ${dir_out}/xxxx.fasta

# Use FIMO to find individual occurences of motif
mkdir ${dir_out}/results/fimo
fimo -oc ${dir_out}/results/fimo \
${dir_motifs}/raw_data/meme_motif_curated/meme_motif_curated.meme \
${workdir}/results/fasta/CD83.fasta