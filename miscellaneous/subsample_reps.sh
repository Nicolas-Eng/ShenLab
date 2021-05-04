#!/bin/bash
#$ -l h_rt=20:0:0
#$ -l mem_free=20G
#$ -S /bin/bash
#$ -cwd
#$ -j y


conda activate seqtk 
ifile=$1
seed=$2
percentage=$3
ofile=$4

seqtk sample -s$seed $ifile $percentage > $ofile
echo $command
$command
gzip $ofile

#Subsample IJ241 20% 40% 60% 80%, matched seeds for PE



#Subsample IJ242 20% 40% 60% 80%, matched seeds for PE