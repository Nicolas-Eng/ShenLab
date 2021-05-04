library(ArchR);
set.seed(1);
addArchRThreads(threads = 16) 
inputFiles<-'/shen/shenlabstore3/neng/20200509_scATAC_sneddon/Epcam-12wpc-scATAC/EpcamFragments/fragments.tsv.gz'
name<-'Epcam-12wpc-scATAC'
addArchRGenome("hg38")
ArrowFiles<-'/shen/shenlabstore3/neng/20200509_scATAC_sneddon/Epcam.arrow'
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Epcam-12wpc-scATAC",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)
getAvailableMatrices(proj)
proj <- filterDoublets(ArchRProj = proj)

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", 
                        name = "IterativeLSI", force=TRUE, binarize = TRUE, 
                        filterQuantile=0.97, dimsToUse=1:20)

proj <- addClusters(input = proj, reducedDims = "IterativeLSI", force=TRUE)
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", force=TRUE)
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h") 


ENDOonly<-proj[which(proj$Clusters=='C5' | proj$Clusters=='C6' | proj$Clusters=='C7' | proj$Clusters=='C8' | proj$Clusters=='C9')]

saveArchRProject(ArchRProj = ENDOonly, outputDirectory = "ENDOonly", load = FALSE)
