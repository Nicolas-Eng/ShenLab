#!/bin/bash
#$ -l h_rt=336:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y

filename=$1

bamCoverage --bam $1.bam --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --ignoreForNormalization chrX -o $1.bw


