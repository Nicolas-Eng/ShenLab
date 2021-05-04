#!/bin/bash
#$ -l h_rt=10:0:0
#$ -l mem_free=50G
#$ -S /bin/bash
#$ -cwd
#$ -j y


cd /shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/merged

#Merged Replicate bigWig Files
bigWigMerge /shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/Human1_RNA_Mock_24h_merged.plus.bw \
	/shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/Human3_RNA_Mock_24h_merged.plus.bw \
	macrophage.plus.bedGraph

bigWigMerge /shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/Human1_RNA_Mock_24h_merged.minus.bw \
	/shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/Human3_RNA_Mock_24h_merged.minus.bw \
	macrophage.minus.bedGraph

#Sort Bedgraph file
LC_COLLATE=C sort -k1,1 -k2,2n macrophage.plus.bedGraph > macrophage.plus.sorted.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n macrophage.minus.bedGraph > macrophage.minus.sorted.bedGraph

#Bedgraph to BigWig
bedGraphToBigWig macrophage.plus.sorted.bedGraph \
	/shen/shenlabstore3/neng/reference_genome/hg19/rsem_star/GRCh37.chrom.sizes \
	macrophage.plus.transcription.bw

bedGraphToBigWig macrophage.minus.sorted.bedGraph \
	/shen/shenlabstore3/neng/reference_genome/hg19/rsem_star/GRCh37.chrom.sizes \
	macrophage.minus.transcription.bw
 
