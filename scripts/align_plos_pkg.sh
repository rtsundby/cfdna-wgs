#!/usr/bin/env bash
input_ref=$1
input_r1=$2
input_r2=$3
threads=$4
output_sort=$5

module load bwa/0.7.17
module use --append ~/modulefiles
module load samtools/1.7
bwa mem -M -t $threads \
    $input_ref \
    $input_r1 \
    $input_r2 |
    samtools view -@ $threads -Sb - -o - |
    samtools sort -@ $threads - -o $output_sort
samtools index -@ threads $output_sort

#input_ref="/data/sundbyrt/reference/hg19_ref_genome/GCF_000001405.25_GRCh37.p13_genomic.fna.gz"
#input_r1="/data/ShernLiquidBx/cfdna-wgs-data-hg19/lib452_processed_R1.fastq.gz"
#input_r2="/data/ShernLiquidBx/cfdna-wgs-data-hg19/lib452_processed_R2.fastq.gz"
#output_sort="/data/ShernLiquidBx/cfdna-wgs-data-hg19/cfdna-wgs-bams/lib452_raw.bam",
#threads="10"

