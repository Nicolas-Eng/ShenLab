#/bin/bash
#$ -l h_rt=30:0:0
#$ -l mem_free=30G
#$ -S /bin/bash
#$ -cwd
#$ -j y

date
hostname
HOME=/shen/shenlabstore3/neng/20200106_H3K27ac_C+R/analysis
INDEX=$HOME/index.txt
SAMPLE=$1
GENOME=hg38
#GENOME=hg19
#REFERENCE=/netapp/home/neng/reference_genomes/
LANE=/shen/shenlabstore3/neng/20200106_H3K27ac_C+R
conda activate cutnrun

echo $LANE/${SAMPLE}*_R1*.f*.gz

#fastp -t 100 -T 100 -i $LANE/*${SAMPLE}*_R1.fq.gz -I $LANE/*${SAMPLE}*_R2.fq.gz -o $HOME/${SAMPLE}.fastp.50bp.R1.fq.gz -O $HOME/${SAMPLE}.fastp.50bp.R2.fq.gz -p -h $HOME/"${SAMPLE}.fastp.html" -j $HOME/"${SAMPLE}.fastp.json"

# Parallelize samples based on input file in run folder.
echo Processing $SAMPLE now...
input1=$HOME/$SAMPLE.fastp.50bp.R1.fq.gz
input2=$HOME/$SAMPLE.fastp.50bp.R2.fq.gz
echo input_one=${input1}
echo input_two=${input2}

mkdir $HOME/${SAMPLE}_analysis_${GENOME}_bowtie2_10_700
cd $HOME/${SAMPLE}_analysis_${GENOME}_bowtie2_10_700

samtools view -h ${SAMPLE}.nuc.nodups.bam | awk -f /wynton/home/shen/neng/scripts/filter.120.awk | samtools view -Sb - > ${SAMPLE}.120bp.bam
samtools view -h ${SAMPLE}.nuc.nodups.bam | awk -f /wynton/home/shen/neng/scripts/filter.150.awk | samtools view -Sb - > ${SAMPLE}.150bp.bam

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

