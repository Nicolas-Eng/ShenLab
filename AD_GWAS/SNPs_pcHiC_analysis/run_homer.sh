#!/bin/bash
#$ -l h_rt=20:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y

conda activate HOMER

cd /shen/shenlabstore3/neng/20200309_MGEN/SNPs_pcHiC_analysis/scripts/final_MGcomp

mkdir MG_HOMERresults
mkdir MG_InVitro_HOMERresults
mkdir MG_ExVivo_HOMERresults

#findMotifsGenome.pl MG.distalATAC.500bp.srt.merged.bed hg19 MG_HOMERresults -size 300
#findMotifsGenome.pl MG_InVitro.distalATAC.500bp.srt.merged.bed hg19 MG_InVitro_HOMERresults -size 300
findMotifsGenome.pl MG_ExVivo.distalATAC.500bp.srt.merged.bed hg19 MG_ExVivo_HOMERresults -size 300
