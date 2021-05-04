#$ -l h_rt=336:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -t 1-100

cd /shen/shenlabstore3/neng/20200509_scATAC_sneddon/SNP_analysis
conda activate EnhancerSwitch

celltype=$1
#inpeaks=12wpc_endocrine_peaks/${celltype}/${celltype}_Peaks_0.01.bed
inpeaks=12wpc_endocrine_peaks/${celltype}/${celltype}_Peaks_0.1.bed

setnumber=$(($SGE_TASK_ID+1))
ctrlfile=12wpc_endocrine_control/${celltype}_Peaks_ctrl_s${setnumber}.bed


bedtools shuffle -i ${inpeaks} -seed ${setnumber} -g hg38.chrom.sizes > ${ctrlfile}.tmp
cut -f1-3 ${ctrlfile}.tmp > ${ctrlfile}
rm ${ctrlfile}.tmp