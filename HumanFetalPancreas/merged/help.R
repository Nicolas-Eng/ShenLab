library(ArchR);
set.seed(1);

setwd('/shen/shenlabstore3/neng/20200509_scATAC_sneddon/merged')
addArchRThreads(threads = 16) 
addArchRGenome("hg38")


outdir<-'/shen/shenlabstore3/neng/20200509_scATAC_sneddon/merged/Merged_16wpc'
arrowFiles<-c('Ep-Nu_2.arrow','Total-Nu_1.arrow','Total-Nu_2.arrow') 
Merge_16wpc <- ArchRProject(ArrowFiles= arrowFiles, outputDirectory=outdir, copyArrows=TRUE)