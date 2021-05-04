#!/bin/bash
#$ -l h_rt=100:0:0
#$ -l mem_free=80G
#$ -S /bin/bash
#$ -cwd
#$ -j y
cd /shen/shenlabstore3/neng/20200509_scATAC_sneddon
conda activate ArchR.3.6.1
Rscript Epcam_pipeline.R