#!/bin/bash
#$ -l h_rt=10:0:0
#$ -l mem_free=50G
#$ -S /bin/bash
#$ -cwd
#$ -j y


cd /shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/merged

#Merged Replicate bigWig Files
bigWigMerge /shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/LMA03_trim_merged.plus.bw \
	/shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/LMA04_trim_merged.plus.bw \
	MG_day7.plus.bedGraph

bigWigMerge /shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/LMA03_trim_merged.minus.bw \
	/shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/LMA04_trim_merged.minus.bw \
	MG_day7.minus.bedGraph

#Sort Bedgraph file
LC_COLLATE=C sort -k1,1 -k2,2n MG_day7.plus.bedGraph > MG_day7.plus.sorted.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n MG_day7.minus.bedGraph > MG_day7.minus.sorted.bedGraph

#Bedgraph to BigWig
bedGraphToBigWig MG_day7.plus.sorted.bedGraph \
	/shen/shenlabstore3/neng/reference_genome/hg19/rsem_star/GRCh37.chrom.sizes \
	MG_day7.plus.transcription.bw

bedGraphToBigWig MG_day7.minus.sorted.bedGraph \
	/shen/shenlabstore3/neng/reference_genome/hg19/rsem_star/GRCh37.chrom.sizes \
	MG_day7.minus.transcription.bw
 
