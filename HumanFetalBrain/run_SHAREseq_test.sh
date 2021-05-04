#!/bin/bash
#$ -l h_rt=336:0:0
#$ -l mem_free=100G
#$ -S /bin/bash
#$ -cwd
#$ -j y

cd /shen/shenlabstore3/neng/20201116_SHAREseq

conda activate SHAREseq

sh SHAREseq_test_RNA.sh