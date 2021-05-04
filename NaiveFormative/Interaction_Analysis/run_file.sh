#!/bin/bash
#$ -l h_rt=30:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y 


conda activate EnhancerSwitch
cd /shen/shenlabstore3/neng/20200528_EnhancerSwitch/Interaction_Analysis

Rscript gene_interactions.R