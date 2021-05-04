#!/bin/bash
#$ -l h_rt=336:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -t 1-2

i=$(($SGE_TASK_ID - 1))
Prefix_list=("naive" "formative")

file_prefix=${Prefix_list[i]}

cd /shen/shenlabstore3/neng/20200617_NaiveForm/RNA-seq/output/visualization/merged

#Merged Replicate bigWig Files
bigWigMerge ../${file_prefix}_RNA_rep1.plus.bw \
	../${file_prefix}_RNA_rep2.plus.bw \
	../${file_prefix}_RNA_rep3.plus.bw \
	${file_prefix}.plus.bedGraph

bigWigMerge ../${file_prefix}_RNA_rep1.minus.bw \
	../${file_prefix}_RNA_rep2.minus.bw \
	../${file_prefix}_RNA_rep3.minus.bw \
	${file_prefix}.minus.bedGraph

#Sort Bedgraph file
LC_COLLATE=C sort -k1,1 -k2,2n ${file_prefix}.plus.bedGraph > ${file_prefix}.plus.sorted.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n ${file_prefix}.minus.bedGraph > ${file_prefix}.minus.sorted.bedGraph

#Bedgraph to BigWig
bedGraphToBigWig ${file_prefix}.plus.sorted.bedGraph \
	/shen/shenlabstore3/neng/reference_genome/mm10/rsem_star_new/GRCm38.p6.sizes.genome \
	${file_prefix}.plus.transcription.bw

bedGraphToBigWig ${file_prefix}.minus.sorted.bedGraph \
	/shen/shenlabstore3/neng/reference_genome/mm10/rsem_star_new/GRCm38.p6.sizes.genome \
	${file_prefix}.minus.transcription.bw
 
