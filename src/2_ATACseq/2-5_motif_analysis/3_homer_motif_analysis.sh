#!/bin/sh

# specify input and output directories
dir_bed="directory containing input bed file/xxxx.bed"
dir_out="output directory"
dir_motifs="directory containing motifs obtained from Michida et al., 2020"

# Specify fasta (fa) directory 
dir_fasta="directory containing fasta file/xxxx.fa"

# Find motifs using Homer using motifs files from Michida et al., 2020
for i in PU1 NFkB
do
	for k in SE_gain SE_unchanged SE_lost TE_gain TE_unchanged TE_lost
	do
        findMotifsGenome.pl ${dir_bed}/${k}.bed \
        ${dir_fasta}/GRCg6a_all.fa \
        ./${dir_out} -size given -p 64 -find ${dir_motifs}/${i}_homer.motif > ${dir_out}/${i}_${k}.txt
	done
done
