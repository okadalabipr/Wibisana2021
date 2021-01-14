#!/bin/sh

# Merge bed file of stimulated and unstimulated cells

# Specify fasta (fa) and annotation (GTF) directories 
dir_fasta="directory containing fasta file/xxxx.fa"
dir_gtf="directory containing gtf file/xxxx.gtf"

# specify input and output directories
dir_out="output directory"
dir_stim="homer tag directory of stimulated cells (all peaks)"
dir_unstim="homer tag directory of stimulated cells (all peaks)"
dir_stim_5k="homer tag directory of stimulated cells (5kb stitched)"
dir_unstim_5k="homer tag directory of stimulated cells (5kb stitched)"

mkdir ${dir_out}/combined_bed
mkdir ${dir_out}/combined_bed_annotated

# combine all data to get all peaks without stitching
cat ${dir_stim}/superEnhancers.bed \
${dir_stim}/typicalEnhancers.bed \
${dir_unstim}/superEnhancers.bed \
${dir_unstim}/typicalEnhancers.bed | \
sed '/^#/d' | \
sort -k1,1 -k2,2n | \
mergeBed > ${dir_out}/combined_bed/all_peaks_unstitched.bed

# merge super enhancers between unstim and stim
cat ${dir_stim_5k}superEnhancers.bed \
${dir_unstim_5k}/superEnhancers.bed | \
sed '/^#/d' | \
sort -k1,1 -k2,2n | \
mergeBed > ${dir_out}/combined_bed/SE_5000_combined.bed

# merge typical enhancers between unstim and stim
cat ${dir_stim_5k}/typicalEnhancers.bed \
${dir_unstim_5k}/typicalEnhancers.bed | \
sed '/^#/d' | \
sort -k1,1 -k2,2n | \
mergeBed > ${dir_out}/combined_bed/TE_5000_combined.bed

# # get constituent peaks at super enhancers 
# intersectBed -a ${workdir}/results/combined_bed/all_peaks_unstitched.bed \
# -b ${workdir}/results/combined_bed/SE_5000_combined.bed > ${workdir}/results/combined_bed/SE_constituent.bed

# # get constituent peaks at typical enhancers
# intersectBed -a ${workdir}/results/combined_bed/all_peaks_unstitched.bed \
# -b ${workdir}/results/combined_bed/TE_5000_combined.bed > ${workdir}/results/combined_bed/TE_constituent.bed

# Annotate all peaks
for i in all_peaks_unstitched SE_5000_combined TE_5000_combined
do
	annotatePeaks.pl ${dir_out}/combined_bed/${i}.bed \
	${dir_fasta} \
	-gtf ${dir_gtf} \
	-norm 1000000 \
	-size given > ${dir_out}/combined_bed_annotated/${i}_annotated.txt
done

# Annotate stitched peaks with intensity (perform for both tag directories, stimulated and unstimulated cells)
for i in SE_5000_combined TE_5000_combined
do
	annotatePeaks.pl ${workdir}/results/combined_bed/${i}.bed \
	${dir_fasta} \
	-gtf ${dir_gtf} \
	-norm 1000000 \
	-size given \
	-d ${dir_out}/tag_directory ${dir_out}/tag_directory > ${dir_out}/combined_bed_annotated/${i}_annotated_intensity.txt

done



