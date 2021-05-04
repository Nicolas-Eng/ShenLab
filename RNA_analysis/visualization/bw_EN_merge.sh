#!/bin/bash
#$ -l h_rt=10:0:0
#$ -l mem_free=50G
#$ -S /bin/bash
#$ -cwd
#$ -j y


cd /shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/merged

#Merged Replicate bigWig Files
bigWigMerge /shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/LM007_trim_merged.plus.bw \
	/shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/LM008_trim_merged.plus.bw \
	EN.plus.bedGraph

bigWigMerge /shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/LM007_trim_merged.minus.bw \
	/shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/LM008_trim_merged.minus.bw \
	EN.minus.bedGraph

#Sort Bedgraph file
LC_COLLATE=C sort -k1,1 -k2,2n EN.plus.bedGraph > EN.plus.sorted.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n EN.minus.bedGraph > EN.minus.sorted.bedGraph

#Bedgraph to BigWig
bedGraphToBigWig EN.plus.sorted.bedGraph \
	/shen/shenlabstore3/neng/reference_genome/hg19/rsem_star/GRCh37.chrom.sizes \
	EN.plus.transcription.bw

bedGraphToBigWig EN.minus.sorted.bedGraph \
	/shen/shenlabstore3/neng/reference_genome/hg19/rsem_star/GRCh37.chrom.sizes \
	EN.minus.transcription.bw
 
