library(ArchR)
inputFiles<-'/shen/shenlabstore3/neng/20200509_scATAC_sneddon/EpcamFragments/outs/fragments.tsv.gz'
name<-'Epcam'
addArchRGenome("hg38")
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = name,
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
addArchRGenome("hg38")
addArchRThreads(threads = 16) 

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)
projEpcam1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "EpcamTutorial",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
projEpcam1
saveArchRProject(ArchRProj = projEpcam1, outputDirectory = "Save-ProjEpcam1", load = FALSE)

projEpcam1<-loadArchRProject(path='/shen/shenlabstore3/neng/20200509_scATAC_sneddon/Save-ProjEpcam1',force=TRUE, showLogo=TRUE)
projEpcam2<-projEpcam1
projEpcam2 <- addIterativeLSI(
    ArchRProj = projEpcam2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.7), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:20
)

projEpcam2 <- addClusters(
    input = projEpcam2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)
projEpcam2 <- addUMAP(
    ArchRProj = projEpcam2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)
saveArchRProject(ArchRProj = projEpcam2, outputDirectory = "Save-ProjEpcam2", load = FALSE)

projEpcam2<-loadArchRProject(path='/shen/shenlabstore3/neng/20200509_scATAC_sneddon/Save-ProjEpcam2',force=TRUE, showLogo=TRUE)
p1 <- plotEmbedding(ArchRProj = projEpcam2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projEpcam2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

ggAlignPlots(p1, p2, type = "h")
markersGS <- getMarkerFeatures(
    ArchRProj = projEpcam2, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
save(markersGS,file='/shen/shenlabstore3/neng/20200509_scATAC_sneddon/Save-ProjEpcam2/markerGS.Rdata')
load('/shen/shenlabstore3/neng/20200509_scATAC_sneddon/Save-ProjEpcam2/markerGS.Rdata')
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.
01 & Log2FC >= 1.25")

markerGenes <-  c('GAPDH','ACTB','INS','INS-IGF2','SCG5','NKX6-1',
  'PDX1','GCG', 'IRX2', 'IRX1', 'ARX',
  'MAFB', 'TTR','SST','PPY','GHRL',
  'NEUROG3','FEV','SUSD2','CPA1',
  'PTF1A', 'CPA2','REG1A','KRT19',
  'CFTR','SOX9','SPP1','VIM',
  'COL1A1','COL3A1','PTPRC','RAC2',
  'SOX10','WT1',
  'CAV1','PECAM1','EPCAM')


heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

projEpcam2 <- addImputeWeights(projEpcam2)
p <- plotEmbedding(
    ArchRProj = projEpcam2, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projEpcam2)
)
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    tEpcam_ArchR(baseSize = 6.5) +
    tEpcam(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    tEpcam(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
    ArchRProj = projEpcam2, 
    addDOC = FALSE, width = 5, height = 5)


p <- plotBrowserTrack(
    ArchRProj = projEpcam2, 
    groupBy = "Clusters", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$CD34)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes.pdf", 
    ArchRProj = projEpcam2, 
    addDOC = FALSE, width = 5, height = 5)

