#!/bin/bash
#$ -l h_rt=336:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -t 1-8

i=$(($SGE_TASK_ID - 1))
#Prefix_list=("epiLCs_RNA_rep1" "epiLCs_RNA_rep2" "epiLCs_RNA_rep3" \
#	"mESCs_RNA_rep1" "mESCs_RNA_rep2" "mESCs_RNA_rep3")

#Prefix_list=("ESC_rep1" "ESC_rep2" "EpiLC_minusActivin" "EpiLC_plusActivin")

#file_prefix=${Prefix_list[i]}

cd /shen/shenlabstore3/neng/20200528_EnhancerSwitch/Interaction_Analysis
atacbed_list=("mESCs.optimal.atac.bed" "mESCs.optimal.atac.bed" \
	"mESCs.optimal.atac.bed" "mESCs.optimal.atac.bed" \
	"epiLCs.optimal.atac.bed" "epiLCs.optimal.atac.bed" \
	"epiLCs.optimal.atac.bed" "epiLCs.optimal.atac.bed")


distalbed_list=("mESCs.decrease.decrease.distal.bed" "mESCs.increase.decrease.distal.bed" \
	"mESCs.decrease.increase.distal.bed" "mESCs.increase.increase.distal.bed" \
	"epiLCs.decrease.decrease.distal.bed" "epiLCs.increase.decrease.distal.bed" \
	"epiLCs.decrease.increase.distal.bed" "epiLCs.increase.increase.distal.bed")

outbed_list=("mESCs.decrease.decrease.x.atac.bed" "mESCs.increase.decrease.x.atac.bed" \
	"mESCs.decrease.increase.x.atac.bed" "mESCs.increase.increase.x.atac.bed" \
	"epiLCs.decrease.decrease.x.atac.bed" "epiLCs.increase.decrease.x.atac.bed" \
	"epiLCs.decrease.increase.x.atac.bed" "epiLCs.increase.increase.x.atac.bed")

atacbed=${atacbed_list[i]}
distalbed=${distalbed_list[i]}
outbed=${outbed_list[i]}


#bedtools intersect -a peaks/atac/${atacbed} -b downstream_new/${distalbed} \
#	| sort -k1,1 -k2,2n - | bedtools merge -i - > downstream_new/output/${outbed}


bedtools intersect -a peaks/atac/${atacbed} -b downstream_new_0.000001/${distalbed} \
	| sort -k1,1 -k2,2n - | bedtools merge -i - > downstream_new_0.000001/output/${outbed}