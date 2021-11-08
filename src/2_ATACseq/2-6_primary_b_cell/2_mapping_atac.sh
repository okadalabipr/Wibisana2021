#!/bin/sh

# edit samplesheet.csv as described in https://nf-co.re/atacseq/1.2.1/usage

# run pipeline to produce bw file
nextflow run nf-core/atacseq -r 1.2.1 \
    --input samplesheet.csv \
    --genome GRCm38 \
    -profile singularity \
    --save_reference \
    --max_cpus 64 \
    --max_memory 256.GB


