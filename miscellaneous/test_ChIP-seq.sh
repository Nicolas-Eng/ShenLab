
INPUT=/netapp/home/neng/scripts/test_ChIP-seq.json
metadata=/netapp/home/neng/scripts/05202019metadata.json
java -jar -Dconfig.file=/netapp/home/neng/dependencies/chip-seq-pipeline2/backends/backend.conf /netapp/home/neng/dependencies/chip-seq-pipeline2/cromwell-34.jar run /netapp/home/neng/dependencies/chip-seq-pipeline2/chip.wdl -i ${INPUT} -m ${PIPELINE_METADATA}