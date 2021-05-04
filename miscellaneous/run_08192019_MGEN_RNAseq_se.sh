#!/bin/bash
#$ -l h_rt=336:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y


file_prefix=$1


conda activate rna-seq-pipeline

echo 'running trimgalore'
trim_galore -q 20 --phred33 --gzip --length 20 --stringency 3 \
	--trim-n --output_dir /shen/shenlabstore3/neng/08192019_MGEN/RNA_analysis/output/trimmed \
	/shen/shenlabstore3/neng/08192019_MGEN/RNA_analysis/data/${file_prefix}_R1.fastq.gz

echo 'running star...'
STAR --genomeDir /shen/shenlabstore3/neng/reference_genomes/rsem_star/hg19 \
	--readFilesCommand zcat --readFilesIn \
	/shen/shenlabstore3/neng/08192019_MGEN/RNA_analysis/output/trimmed/${file_prefix}_R1_trimmed.fq.gz \
	--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within \
	--outFilterMultimapNmax 20 --quantMode TranscriptomeSAM --runThreadN 8 --outFileNamePrefix \
	/shen/shenlabstore3/neng/08192019_MGEN/RNA_analysis/output/${file_prefix}. --outFilterType BySJout \
	--outSAMattributes NH HI AS NM MD --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
	--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 --sjdbScore 1
	
echo 'running rsem...'
rsem-calculate-expression --bam --no-bam-output -p 8 --forward-prob 0 --estimate-rspd \
	/shen/shenlabstore3/neng/08192019_MGEN/RNA_analysis/output/${file_prefix}.Aligned.toTranscriptome.out.bam \
	/shen/shenlabstore3/neng/reference_genomes/rsem_star/hg19/hg19_gencode \
	/shen/shenlabstore3/neng/08192019_MGEN/RNA_analysis/output/${file_prefix} >& \
	/shen/shenlabstore3/neng/08192019_MGEN/RNA_analysis/output/${file_prefix}.rsem.log

echo 'preparing visualization files'

cd /shen/shenlabstore3/neng/08192019_MGEN/RNA_analysis/output
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
/shen/shenlabstore3/neng/reference_genomes/rsem_star/hg19/GRCh37.chrom.sizes \
${file_prefix}.plus.bw
bedGraphToBigWig ${file_prefix}.minus.sorted.bedGraph \
/shen/shenlabstore3/neng/reference_genomes/rsem_star/hg19/GRCh37.chrom.sizes \
${file_prefix}.minus.bw