library(ArchR);
set.seed(1);

setwd('/shen/shenlabstore3/neng/20200509_scATAC_sneddon/merged')
addArchRThreads(threads = 16) 
addArchRGenome("hg38")
outdir<-'/shen/shenlabstore3/neng/20200509_scATAC_sneddon/merged/revision3'
inputFiles<-c('/shen/shenlabstore3/neng/20200509_scATAC_sneddon/merged/fragments/Epcam-12wpc.fragments.tsv.gz',
    '/shen/shenlabstore3/neng/20200509_scATAC_sneddon/merged/fragments/Ep-Nu_2.fragments.tsv.gz',
    '/shen/shenlabstore3/neng/20200509_scATAC_sneddon/merged/fragments/Total-Nu_1.fragments.tsv.gz',
    '/shen/shenlabstore3/neng/20200509_scATAC_sneddon/merged/fragments/Total-Nu_2.fragments.tsv.gz')

names(inputFiles)<-c('Epcam-12wpc','Ep-Nu_2','Total-Nu_1','Total-Nu_2')
proj <- loadArchRProject(path = outdir)

markerGenes <-  c('GAPDH','ACTB','INS','INS-IGF2','SCG5','NKX6-1',
  'PDX1','GCG', 'IRX2', 'IRX1', 'ARX',
  'MAFB', 'TTR','SST','PPY','GHRL',
  'NEUROG3','FEV','SUSD2','CPA1',
  'PTF1A', 'CPA2','REG1A','KRT19',
  'CFTR','SOX9','SPP1','VIM',
  'COL1A1','COL3A1','PTPRC','RAC2',
  'SOX10','WT1',
  'CAV1','PECAM1','EPCAM')

proj <- addImputeWeights(proj)
p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)

plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)
