#!/bin/sh

annotatePeaks.pl conns.bed \
    GRCg6a_all.fa \
    -gtf Gallus_gallus.GRCg6a.96.gtf \
    -size given > conns_annotated.txt
