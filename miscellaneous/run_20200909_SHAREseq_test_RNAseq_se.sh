#!/bin/bash
#$ -l h_rt=336:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -t 1

#i=$(($SGE_TASK_ID - 1))
#Prefix_list=("epiLCs_RNA_rep1" "epiLCs_RNA_rep2" "epiLCs_RNA_rep3" \
#	"mESCs_RNA_rep1" "mESCs_RNA_rep2" "mESCs_RNA_rep3")

#Prefix_list=("ESC_rep1" "ESC_rep2" "EpiLC_minusActivin" "EpiLC_plusActivin")

#file_prefix=${Prefix_list[i]}


conda activate rna-seq-pipeline

echo 'running trimgalore'
#trim_galore -q 20 --phred33 --gzip --length 20 --stringency 3 \
#	--trim-n --output_dir /shen/shenlabstore3/neng/20200106_mESCs_epiLCs/RNA-seq/output/trimmed \
#	/shen/shenlabstore3/neng/20200106_mESCs_epiLCs/RNA-seq/data/${file_prefix}_R1.fastq

echo 'running star...'
STAR --genomeDir /shen/shenlabstore3/neng/reference_genome/mm10/rsem_star_new \
	--readFilesCommand zcat --readFilesIn \
	/shen/shenlabstore3/neng/20200909_SHAREseq_test/polyT/TCGGACGA+TAGATCGC_polyT_R1.fastq.gz \
	--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within \
	--outFilterMultimapNmax 20 --quantMode TranscriptomeSAM --runThreadN 8 --outFileNamePrefix \
	/shen/shenlabstore3/neng/20200909_SHAREseq_test/STAR_out. --outFilterType BySJout \
	--outSAMattributes NH HI AS NM MD --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
	--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 --sjdbScore 1