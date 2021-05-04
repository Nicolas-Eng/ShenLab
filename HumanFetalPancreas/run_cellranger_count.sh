#!/bin/bash
#$ -l h_rt=80:0:0
#$ -l mem_free=100G
#$ -S /bin/bash
#$ -cwd
#$ -j y

cd /shen/shenlabstore3/neng/20200509_scATAC_sneddon

conda activate cellranger-atac
cellranger-atac count --id=Epcam-12wpc-scATAC --reference=/shen/shenlabstore3/neng/dependencies/cellranger-atac-1.2.0/refdata-cellranger-atac-GRCh38-1.2.0 \
	--fastqs=/shen/shenlabstore3/neng/20200509_scATAC_sneddon/raw_data \
	--sample=Epcam-12wpc-scATAC \
	--localcores=8 \
	--localmem=64
