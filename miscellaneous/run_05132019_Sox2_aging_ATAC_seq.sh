#!/bin/bash
#$ -l h_rt=336:0:0
#$ -l mem_free=100G
#$ -S /bin/bash
#$ -cwd
#$ -j y

conda activate encode-atac-seq-pipeline
INPUT=examples/local/ENCSR356KRQ_subsampled.json
PIPELINE_METADATA=metadata.json
java -jar -Dconfig.file=backends/backend.conf cromwell-38.jar run atac.wdl -i ${INPUT} -m ${PIPELINE_METADATA}
