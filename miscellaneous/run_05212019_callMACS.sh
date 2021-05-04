#!/bin/bash
#$ -l h_rt=336:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y

conda activate MACS
echo 'Experiment: '$1


echo 'cd '$2

cd $2

mkdir $2/MACS2

cp $1_merged.shrt.vip.bed $2/MACS2

cd $2/MACS2

macs2 callpeak -t $1_merged.shrt.vip.bed -f BED -g mm -n $1_q0.1_test -B -q 0.1 --nomodel --extsize 147
macs2 callpeak -t $1_merged.shrt.vip.bed -f BED -g mm -n $1_q0.01_test -B -q 0.01 --nomodel --extsize 147
macs2 callpeak -t $1_merged.shrt.vip.bed -f BED -g mm -n $1_q0.05_test -B -q 0.05 --nomodel --extsize 147
macs2 callpeak -t $1_merged.shrt.vip.bed -f BED -g mm -n $1_q0.001_test -B -q 0.001 --nomodel --extsize 147


