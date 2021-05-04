#!/bin/bash
#$ -l h_rt=336:0:0
#$ -l mem_free=100G
#$ -S /bin/bash
#$ -cwd
#$ -j y


cd /shen/shenlabstore3/neng/20200309_MGEN/scripts
conda activate pchic-brain

#Rscript pchic_ne_exploratory.R
Rscript pchic_ne_MGcomp.R