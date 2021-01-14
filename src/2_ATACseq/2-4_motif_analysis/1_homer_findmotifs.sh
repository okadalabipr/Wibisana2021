#!/bin/sh

# list directories
dir_raw="directory containing bed files of gained SE and TE" 
dir_out="output directory"

# Specify fasta (fa) and annotation (GTF) directories 
dir_fasta="directory containing fasta file/xxxx.fa"
dir_gtf="directory containing gtf file/xxxx.gtf"

# Specify CPU and memory to use
N_thread="number of threads"
N_mem="GB of memory"

########### find motifs for intersections ##########
# mkdir ${workdir}/results/SE_activated/

# findMotifsGenome.pl \
# ${workdir}/raw_data/SE_activated_homer.bed \
# ${ref_dir}/GRCg6a_all.fa \
# ${workdir}/results/SE_activated/ \
# -size given \
# -p 16

############ find motifs for SE_gain only ########
mkdir ${dir_out}/results/SE_gain/

findMotifsGenome.pl \
${dir_raw}/SE_gain.bed \
${dir_fasta} \
${dir_out}/results/SE_gain/ \
-size given \
-p ${N_thread}

############ find motifs for TE_gain only ########
mkdir ${dir_out}/results/TE_gain/

findMotifsGenome.pl \
${dir_raw}/TE_gain.bed \
${dir_fasta} \
${dir_out}/results/TE_gain/ \
-size given \
-p ${N_thread}
