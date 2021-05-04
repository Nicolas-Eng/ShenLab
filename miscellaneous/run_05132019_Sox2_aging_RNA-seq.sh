#!/bin/bash
#$ -l h_rt=336:0:0
#$ -l mem_free=32G
#$ -S /bin/bash
#$ -cwd
#$ -j y

conda activate salmon 
#salmon quant -i /shen/shenlabstore3/neng/Salmon_Index/gencode.vM22_salmon-0.8.1-0 --libType A \
#	--gcBias --biasSpeedSamp 5 \
#	-1 /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/data/RNA006_12mSox2hi_1.fq.gz \
#	-2 /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/data/RNA006_12mSox2hi_2.fq.gz \
#	 -o /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/output/Salmon_output/RNA006_12mSox2hi

#salmon quant -i /shen/shenlabstore3/neng/Salmon_Index/gencode.vM22_salmon-0.8.1-0 --libType A \
#	--gcBias --biasSpeedSamp 5 \
#	-1 /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/data/RNA003_2mSox2hi_1.fq.gz \
#	-2 /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/data/RNA003_2mSox2hi_2.fq.gz \
#	 -o /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/output/Salmon_output/RNA003_2mSox2hi

#salmon quant -i /shen/shenlabstore3/neng/Salmon_Index/gencode.vM22_salmon-0.8.1-0 --libType A \
#	--gcBias --biasSpeedSamp 5 \
#	-1 /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/data/RNA009_18mSox2hi_1.fq.gz \
#	-2 /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/data/RNA009_18mSox2hi_2.fq.gz \
#	 -o /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/output/Salmon_output/RNA009_18mSox2hi

#salmon quant -i /shen/shenlabstore3/neng/Salmon_Index/gencode.vM22_salmon-0.8.1-0 --libType A \
#	--gcBias --biasSpeedSamp 5 \
#	-1 /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/data/RNA012_24mSox2hi_1.fq.gz \
#	-2 /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/data/RNA012_24mSox2hi_2.fq.gz \
#	 -o /shen/shenlabstore3/neng/05132019_Sox2_aging/RNA/output/Salmon_output/RNA012_24mSox2hi
