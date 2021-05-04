library(DESeq2)
library(tximport)
library(pheatmap)
library(RColorBrewer)

MasterFile <- read.csv("DESeq_allReps_noEN.csv")
dir <- ("/shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/DESeq")
files <- file.path(dir, paste0(MasterFile$ID, ".merged.genes.results"))
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE,countsCol="expected_count")

names.files <- c("GFiMG Day2 rep1","GFiMG Day2 rep2",
	"GFiMG Day7 rep1","GFiMG Day7 rep2",
	"GFiMG Day14 rep1","GFiMG Day14 rep2",
	"GFiMG Day30 rep1","GFiMG Day30 rep2",
	"Ex vivo microglia (Gosselin et al. 2017) rep1","Ex vivo microglia (Gosselin et al. 2017) rep2","Ex vivo microglia (Gosselin et al. 2017) rep3",
	"In vitro Day7 microglia (Gosselin et al. 2017) rep1","In vitro Day7 microglia (Gosselin et al. 2017) rep2",
	"monocyte-derived macrophage (Heinz et al, 2018) rep1","monocyte-derived macrophage (Heinz et al, 2018) rep2",
	"monocyte-derived macrophage (Carlin et al, 2018) rep1","monocyte-derived macrophage (Carlin et al, 2018) rep2", 
	"WTC iPSC rep1", "WTC iPSC rep2")
# rename column names

#

coldata <- data.frame(c("GFiMG Day2","GFiMG Day2",
	"GFiMG Day7","GFiMG Day7",
	"GFiMG Day14","GFiMG Day14",
	"GFiMG Day30","GFiMG Day30",
	"Ex vivo microglia (Gosselin et al. 2017)","Ex vivo microglia (Gosselin et al. 2017)","Ex vivo microglia (Gosselin et al. 2017)",
	"In vitro Day7 microglia (Gosselin et al. 2017)","In vitro Day7 microglia (Gosselin et al. 2017)",
	"monocyte-derived macrophage (Heinz et al, 2018)","monocyte-derived macrophage (Heinz et al, 2018)",
	"monocyte-derived macrophage (Carlin et al, 2018)","monocyte-derived macrophage (Carlin et al, 2018)", 
	"WTC iPSC", "WTC iPSC"))
colnames(coldata)[1] <- "condition"
coldata$sample <- colnames(txi.rsem$counts)
row.names(coldata) <- colnames(txi.rsem$counts)

# gene cannot have 0 length
expDesign <- coldata
txi.rsem$length[txi.rsem$length == 0] <- 1
library(DESeq2)
ddsTxi <- DESeqDataSetFromTximport(txi.rsem, MasterFile, ~condition)


ddsTxi$condition<- as.factor(c("GFiMG Day2 rep1","GFiMG Day2 rep2",
	"GFiMG Day7 rep1","GFiMG Day7 rep2",
	"GFiMG Day14 rep1","GFiMG Day14 rep2",
	"GFiMG Day30 rep1","GFiMG Day30 rep2",
	"Ex vivo microglia (Gosselin et al. 2017) rep1","Ex vivo microglia (Gosselin et al. 2017) rep2","Ex vivo microglia (Gosselin et al. 2017) rep3",
	"In vitro Day7 microglia (Gosselin et al. 2017) rep1","In vitro Day7 microglia (Gosselin et al. 2017) rep2",
	"monocyte-derived macrophage (Heinz et al, 2018) rep1","monocyte-derived macrophage (Heinz et al, 2018) rep2",
	"monocyte-derived macrophage (Carlin et al, 2018) rep1","monocyte-derived macrophage (Carlin et al, 2018) rep2", 
	"WTC iPSC rep1", "WTC iPSC rep2"))

ddsTxi$ID<-as.factor(c("GFiMG Day2","GFiMG Day2",
	"GFiMG Day7","GFiMG Day7",
	"GFiMG Day14","GFiMG Day14",
	"GFiMG Day30","GFiMG Day30",
	"Ex vivo microglia (Gosselin et al. 2017)","Ex vivo microglia (Gosselin et al. 2017)","Ex vivo microglia (Gosselin et al. 2017)",
	"In vitro Day7 microglia (Gosselin et al. 2017)","In vitro Day7 microglia (Gosselin et al. 2017)",
	"monocyte-derived macrophage (Heinz et al, 2018)","monocyte-derived macrophage (Heinz et al, 2018)",
	"monocyte-derived macrophage (Carlin et al, 2018)","monocyte-derived macrophage (Carlin et al, 2018)", 
	"WTC iPSC", "WTC iPSC"))
#PCA analysis

ddsTxi<-DESeq(ddsTxi)
vst <- vst(ddsTxi)
rlog <- rlog(ddsTxi)
png('RNA-seq_pca_rlog.png',res=300,units="in",width=8,height=8)
plotPCA(rlog,intgroup = "ID", ntop = 500)
dev.off()
png('RNA-seq_pca_vst.png',res=300,units="in",width=8,height=8)
plotPCA(vst, intgroup = "ID", ntop = 500)
dev.off()



# Make a heatmap
sampleDists <- dist(t(assay(vst)))
sampleDists
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rlog$condition
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png('RNA-seq_heatmap_rlog.png',res=300,units="in",width=8,height=8)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, ntop=500)
dev.off()


