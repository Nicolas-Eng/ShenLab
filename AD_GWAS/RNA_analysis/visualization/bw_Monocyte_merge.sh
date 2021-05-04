#!/bin/bash
#$ -l h_rt=10:0:0
#$ -l mem_free=50G
#$ -S /bin/bash
#$ -cwd
#$ -j y


cd /shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/merged

#Merged Replicate bigWig Files
bigWigMerge /shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/MDM-RNAseq-Mock-06h-r3_merged.plus.bw \
	/shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/MDM-RNAseq-Mock-06h-r4_trim_merged.plus.bw \
	monocyte.plus.bedGraph

bigWigMerge /shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/MDM-RNAseq-Mock-06h-r3_merged.minus.bw \
	/shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/output/visualization/MDM-RNAseq-Mock-06h-r4_merged.minus.bw \
	monocyte.minus.bedGraph

#Sort Bedgraph file
LC_COLLATE=C sort -k1,1 -k2,2n monocyte.plus.bedGraph > monocyte.plus.sorted.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n monocyte.minus.bedGraph > monocyte.minus.sorted.bedGraph

#Bedgraph to BigWig
bedGraphToBigWig monocyte.plus.sorted.bedGraph \
	/shen/shenlabstore3/neng/reference_genome/hg19/rsem_star/GRCh37.chrom.sizes \
	monocyte.plus.transcription.bw

bedGraphToBigWig monocyte.minus.sorted.bedGraph \
	/shen/shenlabstore3/neng/reference_genome/hg19/rsem_star/GRCh37.chrom.sizes \
	monocyte.minus.transcription.bw
 
