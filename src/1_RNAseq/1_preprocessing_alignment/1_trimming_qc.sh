#!/bin/sh

# Specify input and output directories
dir_fastq="fastqpath/xxxx.fastq"
dir_out="output directory"
dir_genome="directory of genome"

# Specify fasta (fa) and annotation (GTF) directories 
dir_fasta="directory containing fasta file/xxxx.fa"
dir_gtf="directory containing gtf file/xxxx.gtf"

# specify CPU and memory to use
N_thread="number of threads"
N_mem="GB of memory"

# Adapter trimming and quality check
trim_galore --fastqc \
--trim1 ${dir_fastq} \
-o ${dir_out}

# Indexing of genome
STAR --runThreadN 12 --runMode genomeGenerate \
    --genomeDir ${dir_genome} \
    --genomeFastaFiles ${dir_fasta} \
    --sjdbGTFfile ${dir_gtf} \

# Mapping
STAR --runThreadN ${N_thread} --outSAMtype BAM SortedByCoordinate \
	--genomeDir ${dir_genome} \
	--readFilesIn ${dir_fastq} \
	--readFilesCommand zcat \
	--outFileNamePrefix ${dir_out}

# Count aligned data
featureCounts --verbose -T ${N_thread} -t exon -g gene_id \
-a ${dir_gtf} \
-o ${dir_out}/counts.txt ${dir_out}/*.bam

