#!/bin/sh

# download files
for i in DRR202441 DRR202442
do
    fasterq-dump --split-files -e 8 ${i}
done

pigz -5 -p7 *fastq

