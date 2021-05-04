#!/bin/bash
#$ -l h_rt=24:0:0
#$ -l mem_free=15G
#$ -S /bin/bash
#$ -cwd
#$ -j y


cd /shen/shenlabstore3/neng/20200309_MGEN/ATAC_analysis
conda activate encode-atac-seq-pipeline
 
caper run /shen/shenlabstore3/neng/dependencies/atac-seq-pipeline-master/atac.wdl -i /shen/shenlabstore3/neng/20200309_MGEN/ATAC_analysis/Monocyte.json --cromwell /wynton/home/shen/neng/.caper/cromwell_jar/cromwell-47.jar --womtool /wynton/home/shen/neng/.caper/womtool_jar/womtool-47.jar
