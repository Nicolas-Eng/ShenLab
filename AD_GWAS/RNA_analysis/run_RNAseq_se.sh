#!/bin/bash
#$ -l h_rt=336:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -t 1-6

i=$(($SGE_TASK_ID - 1))
Prefix_list=( "naive_RNA_rep1" "naive_RNA_rep2" "naive_RNA_rep3"\
	"formative_RNA_rep1" "formative_RNA_rep2" "formative_RNA_rep3")


file_prefix=${Prefix_list[i]}

cd /shen/shenlabstore3/neng/20200617_NaiveForm/RNA-seq

conda activate rna-seq-pipeline-new

#cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a "A{18}" -a "T{18}" \
#	-m 20 -o processed_data/${file_prefix}_trimmed_R1.fastq data/${file_prefix}_R1.fastq


#STAR --genomeDir /shen/shenlabstore3/neng/reference_genome/mm10/rsem_star_newest --runMode alignReads \
# 	--readFilesIn	/shen/shenlabstore3/neng/20200617_NaiveForm/RNA-seq/processed_data/${file_prefix}_trimmed_R1.fastq \
#	--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within \
#	--runThreadN 8 --outFileNamePrefix /shen/shenlabstore3/neng/20200617_NaiveForm/RNA-seq/output/${file_prefix}. --outFilterType BySJout \
#	--outFilterMultimapNmax 1 --outFilterMismatchNoverReadLmax 0.05 --seedSearchStartLmax 25 --winAnchorMultimapNmax 100


cd /shen/shenlabstore3/neng/20200617_NaiveForm/RNA-seq/output
mkdir visualization
STAR --runMode inputAlignmentsFromBAM \
	--inputBAMfile ${file_prefix}.Aligned.sortedByCoord.out.bam \
	--outWigType bedGraph --outWigStrand Stranded --outWigNorm RPM \
	--outFileNamePrefix visualization/${file_prefix}.
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