#!/bin/bash
#$ -l h_rt=336:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -t 1-6

i=$(($SGE_TASK_ID - 1))
Prefix_list=("epiLCs_RNA_rep1" "epiLCs_RNA_rep2" "epiLCs_RNA_rep3" \
	"mESCs_RNA_rep1" "mESCs_RNA_rep2" "mESCs_RNA_rep3")

#Prefix_list=("ESC_rep1" "ESC_rep2" "EpiLC_minusActivin" "EpiLC_plusActivin")

file_prefix=${Prefix_list[i]}


conda activate rna-seq-pipeline

echo 'running trimgalore'
#trim_galore -q 20 --phred33 --gzip --length 20 --stringency 3 \
#	--trim-n --output_dir /shen/shenlabstore3/neng/20200106_mESCs_epiLCs/RNA-seq/output/trimmed \
#	/shen/shenlabstore3/neng/20200106_mESCs_epiLCs/RNA-seq/data/${file_prefix}_R1.fastq

echo 'running star...'
#STAR --genomeDir /shen/shenlabstore3/neng/reference_genome/mm10/rsem_star_new \
#	--readFilesCommand zcat --readFilesIn \
#	/shen/shenlabstore3/neng/20200106_mESCs_epiLCs/RNA-seq/output/trimmed/${file_prefix}_R1_trimmed.fq.gz \
#	--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within \
#	--outFilterMultimapNmax 20 --quantMode TranscriptomeSAM --runThreadN 8 --outFileNamePrefix \
#	/shen/shenlabstore3/neng/20200106_mESCs_epiLCs/RNA-seq/output/${file_prefix}. --outFilterType BySJout \
#	--outSAMattributes NH HI AS NM MD --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
#	--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 \
#	--alignSJDBoverhangMin 1 --sjdbScore 1
	
echo 'running rsem...'
#rsem-calculate-expression --bam --no-bam-output -p 8 --forward-prob 0 --estimate-rspd \
#	/shen/shenlabstore3/neng/20200106_mESCs_epiLCs/RNA-seq/output/${file_prefix}.Aligned.toTranscriptome.out.bam \
#	/shen/shenlabstore3/neng/reference_genome/mm10/rsem_star_new/mm10_gencode \
#	/shen/shenlabstore3/neng/20200106_mESCs_epiLCs/RNA-seq/output/${file_prefix} >& \
#	/shen/shenlabstore3/neng/20200106_mESCs_epiLCs/RNA-seq/output/${file_prefix}.rsem.log

echo 'preparing visualization files'

cd /shen/shenlabstore3/neng/20200106_mESCs_epiLCs/RNA-seq/output
mkdir visualization
#STAR --runMode inputAlignmentsFromBAM \
#	--inputBAMfile ${file_prefix}.Aligned.sortedByCoord.out.bam \
#	--outWigType bedGraph --outWigStrand Stranded --outWigNorm RPM \
#	--outFileNamePrefix visualization/${file_prefix}.
cd visualization

# Sort bedGraph files.
LC_COLLATE=C sort -k1,1 -k2,2n -o ${file_prefix}.plus.sorted.bedGraph ${file_prefix}.Signal.Unique.str2.out.bg
LC_COLLATE=C sort -k1,1 -k2,2n -o ${file_prefix}.minus.sorted.bedGraph ${file_prefix}.Signal.Unique.str1.out.bg

# Convert bedGraph to bigWig.
bedGraphToBigWig ${file_prefix}.plus.sorted.bedGraph \
/shen/shenlabstore3/neng/reference_genome/mm10/rsem_star_new/GRCm38.p6.sizes.genome \
${file_prefix}.plus.bw
bedGraphToBigWig ${file_prefix}.minus.sorted.bedGraph \
/shen/shenlabstore3/neng/reference_genome/mm10/rsem_star_new/GRCm38.p6.sizes.genome \
${file_prefix}.minus.bw


