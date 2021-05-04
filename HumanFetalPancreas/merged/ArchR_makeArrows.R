library(ArchR);
set.seed(1);

setwd('/shen/shenlabstore3/neng/20200509_scATAC_sneddon/merged')
addArchRThreads(threads = 16) 
addArchRGenome("hg38")
outdir<-'/shen/shenlabstore3/neng/20200509_scATAC_sneddon/merged/revision1'
#inputFiles<-c('/shen/shenlabstore3/neng/20200509_scATAC_sneddon/merged/fragments/Epcam-12wpc.fragments.tsv.gz',
#    '/shen/shenlabstore3/neng/20200509_scATAC_sneddon/merged/fragments/Ep-Nu_2.fragments.tsv.gz',
#    '/shen/shenlabstore3/neng/20200509_scATAC_sneddon/merged/fragments/Total-Nu_1.fragments.tsv.gz',
#    '/shen/shenlabstore3/neng/20200509_scATAC_sneddon/merged/fragments/Total-Nu_2.fragments.tsv.gz')


inputFiles<-c('/shen/shenlabstore3/neng/20200509_scATAC_sneddon/merged/fragments/Ep-Nu_2.fragments.tsv.gz',
    '/shen/shenlabstore3/neng/20200509_scATAC_sneddon/merged/fragments/Total-Nu_1.fragments.tsv.gz',
    '/shen/shenlabstore3/neng/20200509_scATAC_sneddon/merged/fragments/Total-Nu_2.fragments.tsv.gz')

names(inputFiles)<-c('Ep-Nu_2','Total-Nu_1','Total-Nu_2')

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 8, #Dont set this too high because you can always increase later
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  minFrag=10^3.5,
  maxFrag=20000
)



#doubScores <- addDoubletScores(
#  input = ArrowFiles,
#  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#  knnMethod = "LSI", #Refers to the embedding to use for nearest neighbor search.
#  LSIMethod = 1,
#)


proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = outdir,
  copyArrows = FALSE #This is recommended so that you maintain an unaltered copy for later usage.
)

#getAvailableMatrices(proj)
#proj <- filterDoublets(ArchRProj = proj)
#proj <- filterDoublets(proj, filterRatio = 2)
#proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI") 
#proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
#proj <- saveArchRProject(ArchRProj = proj)
