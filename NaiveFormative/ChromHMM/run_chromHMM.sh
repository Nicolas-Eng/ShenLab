#/bin/bash
#$ -l h_rt=72:0:0
#$ -l mem_free=100G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -t 1

i=$(($SGE_TASK_ID - 1))
Prefix_list=("new_noinput" "new" "final")
#Prefix_list=("ESC_rep1" "ESC_rep2" "EpiLC_minusActivin" "EpiLC_plusActivin")

file_prefix=${Prefix_list[i]}


conda activate ChromHMM

cd /shen/shenlabstore3/neng/20200617_NaiveForm/ChromHMM



#java -mx1600M -jar /shen/shenlabstore3/neng/dependencies/ChromHMM/ChromHMM.jar \
#	BinarizeBam /shen/shenlabstore3/neng/20200617_NaiveForm/ChromHMM/mm10.txt \
#	/shen/shenlabstore3/neng/20200617_NaiveForm/ChromHMM/bam \
#	naive_cellmarkfile.txt \
#	naive


#java -mx1600M -jar /shen/shenlabstore3/neng/dependencies/ChromHMM/ChromHMM.jar \
#	BinarizeBam /shen/shenlabstore3/neng/20200617_NaiveForm/ChromHMM/mm10.txt \
#	/shen/shenlabstore3/neng/20200617_NaiveForm/ChromHMM/bam \
#	form_cellmarkfile.txt \
#	form



#java -mx1600M -jar /shen/shenlabstore3/neng/dependencies/ChromHMM/ChromHMM.jar \
#	LearnModel /shen/shenlabstore3/neng/20200617_NaiveForm/ChromHMM/naive \
#	naive/12_states \
#	12 mm10
#### Minimal loglikelihood -7725744.508 #####


#java -mx1600M -jar /shen/shenlabstore3/neng/dependencies/ChromHMM/ChromHMM.jar \
#	LearnModel /shen/shenlabstore3/neng/20200617_NaiveForm/ChromHMM/form \
#	form/12_states \
#	12 mm10
#### Minimal loglikelihood -7241740.700 #####


java -Xmx20g -jar /shen/shenlabstore3/neng/dependencies/ChromHMM/ChromHMM.jar \
   BinarizeBam /shen/shenlabstore3/neng/20200617_NaiveForm/ChromHMM/mm10.txt \
   /shen/shenlabstore3/neng/20200617_NaiveForm/ChromHMM/bam \
   comb_cellmarkfile_${file_prefix}.txt \
   comb_${file_prefix}

java -mx1600M -jar /shen/shenlabstore3/neng/dependencies/ChromHMM/ChromHMM.jar \
	LearnModel -holdcolumnorder -holdroworder /shen/shenlabstore3/neng/20200617_NaiveForm/ChromHMM/comb_${file_prefix} \
	comb_${file_prefix}/21_states_${file_prefix} 21 mm10
