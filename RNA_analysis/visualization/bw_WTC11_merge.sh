#!/bin/bash
#$ -l h_rt=10:0:0
#$ -l mem_free=50G
#$ -S /bin/bash
#$ -cwd
#$ -j y


cd /shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/merged

#Merged Replicate bigWig Files
bigWigMerge /shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/LM124_trim_merged.plus.bw \
	/shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/LM134_trim_merged.plus.bw \
	WTC11.plus.bedGraph

bigWigMerge /shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/LM124_trim_merged.minus.bw \
	/shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/LM134_trim_merged.minus.bw \
	WTC11.minus.bedGraph

#Sort Bedgraph file
LC_COLLATE=C sort -k1,1 -k2,2n WTC11.plus.bedGraph > WTC11.plus.sorted.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n WTC11.minus.bedGraph > WTC11.minus.sorted.bedGraph

#Bedgraph to BigWig
bedGraphToBigWig WTC11.plus.sorted.bedGraph \
	/shen/shenlabstore3/neng/reference_genome/hg19/rsem_star/GRCh37.chrom.sizes \
	WTC11.plus.transcription.bw

bedGraphToBigWig WTC11.minus.sorted.bedGraph \
	/shen/shenlabstore3/neng/reference_genome/hg19/rsem_star/GRCh37.chrom.sizes \
	WTC11.minus.transcription.bw
 
