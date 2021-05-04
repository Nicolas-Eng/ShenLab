#!/bin/bash
#$ -l h_rt=20:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y

conda activate HOMER

cd /shen/shenlabstore3/neng/20200509_scATAC_sneddon/HOMER_analysis

mkdir Beta_HOMERresults
mkdir FEV_HOMERresults
mkdir Common_Progenitor_HOMERresults

#findMotifsGenome.pl Beta.noncoding.bed hg38 Beta_HOMERresults -size 250

#findMotifsGenome.pl FEV.noncoding.bed hg38 FEV_HOMERresults -size 250

findMotifsGenome.pl Common_Progenitor.noncoding.bed hg38 Common_Progenitor_HOMERresults -size 250
