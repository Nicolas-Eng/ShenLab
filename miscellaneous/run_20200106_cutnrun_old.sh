#/bin/bash
#$ -l h_rt=72:0:0
#$ -l mem_free=30G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -t 1-2
date
hostname
HOME=/shen/shenlabstore3/neng/20200106_H3K27ac_C+R/analysis
INDEX=$HOME/index.txt
SAMPLE=$(awk 'FNR == i {print}' i=${SGE_TASK_ID} $INDEX)
GENOME=hg38
#GENOME=hg19
#REFERENCE=/netapp/home/neng/reference_genomes/
LANE=/shen/shenlabstore3/neng/20200106_H3K27ac_C+R
conda activate cutnrun

echo $LANE/${SAMPLE}*_R1*.f*.gz

#fastp -t 100 -T 100 -i $LANE/*${SAMPLE}*_R1*.f*.gz -I $LANE/*${SAMPLE}*_R2*.f*.gz -o $HOME/${SAMPLE}.fastp.50bp.R1.fq.gz -O $HOME/${SAMPLE}.fastp.50bp.R2.fq.gz -p -h $HOME/"${SAMPLE}.fastp.html" -j $HOME/"${SAMPLE}.fastp.json"

# Parallelize samples based on input file in run folder.
echo Processing $SAMPLE now...
input1=$HOME/$SAMPLE*fastp.50bp.R1*.f*
input2=$HOME/$SAMPLE*fastp.50bp.R2*.f*
echo input_one=${input1}
echo input_two=${input2}

mkdir $HOME/${SAMPLE}_analysis_${GENOME}_bowtie2_10_700
cd $HOME/${SAMPLE}_analysis_${GENOME}_bowtie2_10_700

# Perform alignment using bowtie2.
bowtie2  -x /wynton/home/shen/neng/reference_genomes/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 \
	-1 /shen/shenlabstore3/neng/20200106_H3K27ac_C+R/analysis/IJ215.fastp.50bp.R1.fq.gz \
	-2 /shen/shenlabstore3/neng/20200106_H3K27ac_C+R/analysis/IJ215.fastp.50bp.R2.fq.gz \
	| samtools view -bS - > ${SAMPLE}.bam
#bowtie2  --local --very-sensitive-local --no-mixed --no-discordant --phred33 -I 10 -X 700 -x /scrapp2/ijones1/reference_genomes/hg19 -1 ${input1} -2 ${input2} | samtools view -bS - > ${SAMPLE}.bam

# Sort BAM file.
picard SortSam I=${SAMPLE}.bam O=${SAMPLE}.sorted.bam SORT_ORDER=coordinate TMP_DIR=tmp/


# Characterize sorted BAM file.
samtools flagstat ${SAMPLE}.sorted.bam > ${SAMPLE}.sorted.flagstat

# Build BAM index for sorted reads.
picard BuildBamIndex I=${SAMPLE}.sorted.bam TMP_DIR=tmp/

# Process unmapped reads.
samtools view -bhf 4 ${SAMPLE}.sorted.bam > ${SAMPLE}.unmapped.bam
samtools flagstat ${SAMPLE}.unmapped.bam > ${SAMPLE}.unmapped.flagstat
picard BuildBamIndex I=${SAMPLE}.unmapped.bam TMP_DIR=tmp/

# Extract Properly mapped pair
samtools view -bhf 2 ${SAMPLE}.sorted.bam > ${SAMPLE}.proper.pair.bam
samtools flagstat ${SAMPLE}.proper.pair.bam > ${SAMPLE}.proper.pair.flagstat
picard BuildBamIndex I=${SAMPLE}.proper.pair.bam TMP_DIR=tmp/


# Split mitochondrial and nuclear reads.
samtools view -H ${SAMPLE}.proper.pair.bam > ${SAMPLE}.header
samtools view -h ${SAMPLE}.proper.pair.bam | awk '{if($3 == "chrM"){print $0}}' | cat ${SAMPLE}.header - | samtools view -Sb - > ${SAMPLE}.mito.bam
samtools view -h ${SAMPLE}.proper.pair.bam | awk '{if($3 != "chrM"){print $0}}' | samtools view -Sb - > ${SAMPLE}.nuc.bam
rm ${SAMPLE}.header
samtools flagstat ${SAMPLE}.mito.bam > ${SAMPLE}.mito.flagstat
samtools flagstat ${SAMPLE}.nuc.bam > ${SAMPLE}.nuc.flagstat
picard BuildBamIndex I=${SAMPLE}.nuc.bam TMP_DIR=tmp/

# Get unique reads.
picard MarkDuplicates I=${SAMPLE}.nuc.bam O=${SAMPLE}.nuc.nodups.bam REMOVE_DUPLICATES=TRUE METRICS_FILE=${SAMPLE}.markduplicates.txt TMP_DIR=tmp/

samtools view -h ${SAMPLE}.nuc.nodups.bam | awk -f /netapp/home/neng/scripts/filter.120.awk | samtools view -Sb - > ${SAMPLE}.120bp.bam
samtools view -h ${SAMPLE}.nuc.nodups.bam | awk -f /netapp/home/neng/scripts/filter.150.awk | samtools view -Sb - > ${SAMPLE}.150bp.bam

echo "# in ${SAMPLE}.nuc.nodups.bam..."
samtools view  ${SAMPLE}.nuc.nodups.bam | wc -l

echo "# in ${SAMPLE}.120bp.bam..."
samtools view ${SAMPLE}.120bp.bam | wc -l

echo "# in ${SAMPLE}.150bp.bam..."
samtools view ${SAMPLE}.150bp.bam | wc -l

samtools index ${SAMPLE}.nuc.nodups.bam
samtools index ${SAMPLE}.120bp.bam
samtools index ${SAMPLE}.150bp.bam

bamPEFragmentSize \
-hist ${SAMPLE}.fragmentSize.png \
-T "Fragment size of PE C+R data" \
--maxFragmentLength 1000 \
-b ${SAMPLE}.nuc.nodups.bam \
--samplesLabel ${SAMPLE}



# convert bam to bigwig
bamCoverage -b ${SAMPLE}.120bp.bam -of bigwig -o ${SAMPLE}.bowtie2.120.bw
bamCoverage -b ${SAMPLE}.150bp.bam -of bigwig -o ${SAMPLE}.bowtie2.150.bw
bamCoverage -b ${SAMPLE}.nuc.nodups.bam -of bigwig -o ${SAMPLE}.bowtie2.bw

