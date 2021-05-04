#!/bin/bash
#$ -l h_rt=336:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y

filename=$1

samtools sort -o $1.sort.bam $1.bam

samtools index $1.sort.bam

bamCoverage --bam $1.sort.bam --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2150570000 --ignoreForNormalization chrX -o $1.sort.bw


