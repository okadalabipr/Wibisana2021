#!/bin/sh

# perform for both stimulated and unstimulated cells

# specify input and output directories
dir_bam="directory of cellranger-atac output/outs/possorted_bam.bam"
dir_fastq="fastqpath/xxxx.fastq"
dir_out="output directory"
dir_genome="directory of genome"
dir_gtf="directory of gtf file"


mkdir ${dir_out}/results

# make tag directory for normal peak calling
makeTagDirectory ${dir_out}/tag_directory \
${dir_bam} -sspe \
-single \
-tbp 1

# make tag directory for short peak calling
makeTagDirectory ${dir_out}/tag_directory_allpeaks \
${stim_dir}/possorted_bam.bam -sspe \
-single \
-tbp 1

# find peaks using homer
findPeaks ${dir_out}/tag_directory \
-style super -o auto \
-typical ${dir_out}/results/sample_peaks/typicalEnhancers.txt \
-minDist 5000 \
-L 0 -fdr 0.0001

findPeaks ${dir_out}/tag_directory_allpeaks \
-style super -o auto \
-typical ${dir_out}/results/sample_allpeaks/typicalEnhancers.txt \
-minDist 0 \
-L 0 -fdr 0.0001


# annotate peaks
for i in sample_allpeaks sample_peaks
do
	annotatePeaks.pl ${dir_out}/results/${i}/superEnhancers.txt \
    ${dir_genome} \
    -gtf ${dir_gtf} > ${dir_out}/results/${i}/superEnhancer_anno.txt

    annotatePeaks.pl ${dir_out}/results/${i}/typicalEnhancers.txt \
    ${dir_genome} \
    -gtf ${dir_gtf} > ${dir_out}/results/${i}/typicalEnhancer_anno.txt
done

# convert homer position file to conventional bed file 
for i in sample_allpeaks sample_peaks
do
	pos2bed.pl ${dir_out}/results/${i}/superEnhancers.txt > ${dir_out}/results/${i}/superEnhancers.bed
	pos2bed.pl ${dir_out}/results/${i}/typicalEnhancers.txt > ${dir_out}/results/${i}/typicalEnhancers.bed
done