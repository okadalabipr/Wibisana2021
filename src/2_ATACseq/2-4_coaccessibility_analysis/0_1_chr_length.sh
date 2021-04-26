#!/bin/sh

# get chromosome size
samtools faidx genome.fa
cut -f1,2 genome.fa.fai > genomesize.tsv