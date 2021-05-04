#!/bin/bash
#$ -l h_rt=336:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y
conda activate rna-seq-pipeline
echo running trimgalore
trim_galore -q 20 --phred33 --gzip --paired --length 20 --stringency 3 --trim-n --output_dir /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/12m/trimmed /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/12m/RNA006_12mSox2hi_1.fq.gz /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/12m/RNA006_12mSox2hi_2.fq.gz
echo running star...
STAR --genomeDir /shen/shenlabstore3/neng/reference_genomes/rsem_star_mouse/mm10 --readFilesCommand zcat --readFilesIn /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/12m/trimmed/RNA006_12mSox2hi_1_val_1.fq.gz /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/12m/trimmed/RNA006_12mSox2hi_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within --outFilterMultimapNmax 20 --quantMode TranscriptomeSAM --runThreadN 8 --outFileNamePrefix /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/12m/RNA006_12mSox2hi. --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1
echo running rsem...
rsem-calculate-expression --bam --paired-end --no-bam-output -p 8 --forward-prob 0 --estimate-rspd /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/12m/RNA006_12mSox2hi.Aligned.toTranscriptome.out.bam /shen/shenlabstore3/neng/reference_genomes/rsem_star_mouse/mm10/mm10_gencode /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/12m/RNA006_12mSox2hi >& /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/12m/RNA006_12mSox2hi.rsem.log
