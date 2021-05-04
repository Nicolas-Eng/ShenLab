#!/bin/bash
#$ -l h_rt=336:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y 


cd /shen/shenlabstore3/neng/20200617_NaiveForm/RNA-seq

conda activate rna-seq-pipeline-new
outfolder=featureCounts_nostrand

mkdir ${outfolder}

featureCounts -s 0 -t exon -g gene_id -a /shen/shenlabstore3/neng/reference_genome/mm10/rsem_star_newest/gencode.vM25.annotation.gtf \
	-o ${outfolder}/RNA_merged_counts.txt \
	output/naive_RNA_rep1.Aligned.sortedByCoord.out.bam \
	output/naive_RNA_rep2.Aligned.sortedByCoord.out.bam \
	output/naive_RNA_rep3.Aligned.sortedByCoord.out.bam \
	output/formative_RNA_rep1.Aligned.sortedByCoord.out.bam \
	output/formative_RNA_rep2.Aligned.sortedByCoord.out.bam \
	output/formative_RNA_rep3.Aligned.sortedByCoord.out.bam 