#!/bin/sh

# motif analysis for PU.1 and NF-kB

# specify input bed file and motif file
dir_bed="input bed file"
dir_out="output .txt"
dir_genome="genome file"

findMotifsGenome.pl ${dir_bed} \
    ${dir_genome} \
    ./results -size given -p 64 -find ${dir_motif} > ${dir_out}
