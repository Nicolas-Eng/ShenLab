#!/bin/bash
#$ -l h_rt=336:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y

filename=$1

gunzip -f $1.narrowPeak.gz
cut -f1-4 $1.narrowPeak | bedToBam -i - -g /scrapp/neng/genome/genomes/mm10/mm10.chrom.sizes | samtools sort -o $1.sort.bam -

samtools index $1.sort.bam

bamCoverage --bam $1.sort.bam --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2150570000 --ignoreForNormalization chrX -o $1.sort.bw


