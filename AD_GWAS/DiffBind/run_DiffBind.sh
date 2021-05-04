#!/bin/bash
#$ -l h_rt=336:0:0
#$ -l mem_free=100G
#$ -S /bin/bash
#$ -cwd
#$ -j y

conda activate DiffBind

cd /shen/shenlabstore3/neng/20200309_MGEN/DiffBind

Rscript run_DiffBind.R 
