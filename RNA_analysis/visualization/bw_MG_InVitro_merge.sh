#!/bin/bash
#$ -l h_rt=10:0:0
#$ -l mem_free=50G
#$ -S /bin/bash
#$ -cwd
#$ -j y


cd /shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/merged

#Merged Replicate bigWig Files
bigWigMerge /shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/SRR5955149_trim_merged.plus.bw \
	/shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/SRR5955158_trim_merged.plus.bw \
	MG_InVitro.plus.bedGraph

bigWigMerge /shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/SRR5955149_trim_merged.minus.bw \
	/shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/SRR5955158_trim_merged.minus.bw \
	MG_InVitro.minus.bedGraph

#Sort Bedgraph file
LC_COLLATE=C sort -k1,1 -k2,2n MG_InVitro.plus.bedGraph > MG_InVitro.plus.sorted.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n MG_InVitro.minus.bedGraph > MG_InVitro.minus.sorted.bedGraph

#Bedgraph to BigWig
bedGraphToBigWig MG_InVitro.plus.sorted.bedGraph \
	/shen/shenlabstore3/neng/reference_genome/hg19/rsem_star/GRCh37.chrom.sizes \
	MG_InVitro.plus.transcription.bw

bedGraphToBigWig MG_InVitro.minus.sorted.bedGraph \
	/shen/shenlabstore3/neng/reference_genome/hg19/rsem_star/GRCh37.chrom.sizes \
	MG_InVitro.minus.transcription.bw
 
