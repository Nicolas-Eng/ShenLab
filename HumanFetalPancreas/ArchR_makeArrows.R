library(ArchR);
set.seed(1);
addArchRThreads(threads = 16) 
inputFiles<-'/shen/shenlabstore3/neng/20200509_scATAC_sneddon/Epcam-12wpc-scATAC/outs/fragments.tsv.gz'
name<-'Epcam-12wpc-scATAC'
addArchRGenome("hg38")
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = name,
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
