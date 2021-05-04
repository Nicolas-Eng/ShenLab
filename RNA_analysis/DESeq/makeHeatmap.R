# Updated figures for Xioayu 5/4/20

# Read in 800 Microglia Genes of interest ####

load('genesFoundTot.Rdata')
MG_genes <- data.frame(genesFoundTot$ensembl_gene_id)
colnames(MG_genes)[1] <- "gene_id"


# Claire's Day2 #1
MG_Day2_Rep1_genes <- read.table("MG_Day2_Rep1.merged.genes.results", sep="\t", header=TRUE)
MG_Day2_Rep1_genes[,1]<- unlist(strsplit(as.character(MG_Day2_Rep1_genes[, 1]), split="\\."))[(1:length(MG_Day2_Rep1_genes[, 1]))*2-1]
MG_Day2_Rep1_tpm <- MG_Day2_Rep1_genes[,c(1,6)]
colnames(MG_Day2_Rep1_tpm)[2] <- "MG_Day2_Rep1"

# Claire's Day2 #2
MG_Day2_Rep2_genes <- read.table("MG_Day2_Rep2.merged.genes.results", sep="\t", header=TRUE)
MG_Day2_Rep2_genes[,1]<- unlist(strsplit(as.character(MG_Day2_Rep2_genes[, 1]), split="\\."))[(1:length(MG_Day2_Rep2_genes[, 1]))*2-1]
MG_Day2_Rep2_tpm <- MG_Day2_Rep2_genes[,c(1,6)]
colnames(MG_Day2_Rep2_tpm)[2] <- "MG_Day2_Rep2"

# Claire's Day7 #1
MG_Day7_Rep1_genes <- read.table("MG_Day7_Rep1.merged.genes.results", sep="\t", header=TRUE)
MG_Day7_Rep1_genes[,1]<- unlist(strsplit(as.character(MG_Day7_Rep1_genes[, 1]), split="\\."))[(1:length(MG_Day7_Rep1_genes[, 1]))*2-1]
MG_Day7_Rep1_tpm <- MG_Day7_Rep1_genes[,c(1,6)]
colnames(MG_Day7_Rep1_tpm)[2] <- "MG_Day7_Rep1"

# Claire's Day7 #2
MG_Day7_Rep2_genes <- read.table("MG_Day7_Rep2.merged.genes.results", sep="\t", header=TRUE)
MG_Day7_Rep2_genes[,1]<- unlist(strsplit(as.character(MG_Day7_Rep2_genes[, 1]), split="\\."))[(1:length(MG_Day7_Rep2_genes[, 1]))*2-1]
MG_Day7_Rep2_tpm <- MG_Day7_Rep2_genes[,c(1,6)]
colnames(MG_Day7_Rep2_tpm)[2] <- "MG_Day7_Rep2"

# Claire's Day14 #1
MG_Day14_Rep1_genes <- read.table("MG_Day14_Rep1.merged.genes.results", sep="\t", header=TRUE)
MG_Day14_Rep1_genes[,1]<- unlist(strsplit(as.character(MG_Day14_Rep1_genes[, 1]), split="\\."))[(1:length(MG_Day14_Rep1_genes[, 1]))*2-1]
MG_Day14_Rep1_tpm <- MG_Day14_Rep1_genes[,c(1,6)]
colnames(MG_Day14_Rep1_tpm)[2] <- "MG_Day14_Rep1"

# Claire's Day14 #2
MG_Day14_Rep2_genes <- read.table("MG_Day14_Rep2.merged.genes.results", sep="\t", header=TRUE)
MG_Day14_Rep2_genes[,1]<- unlist(strsplit(as.character(MG_Day14_Rep2_genes[, 1]), split="\\."))[(1:length(MG_Day14_Rep2_genes[, 1]))*2-1]
MG_Day14_Rep2_tpm <- MG_Day14_Rep2_genes[,c(1,6)]
colnames(MG_Day14_Rep2_tpm)[2] <- "MG_Day14_Rep2"

# Claire's Day30 #1
MG_Day30_Rep1_genes <- read.table("MG_Day30_Rep1.merged.genes.results", sep="\t", header=TRUE)
MG_Day30_Rep1_genes[,1]<- unlist(strsplit(as.character(MG_Day30_Rep1_genes[, 1]), split="\\."))[(1:length(MG_Day30_Rep1_genes[, 1]))*2-1]
MG_Day30_Rep1_tpm <- MG_Day30_Rep1_genes[,c(1,6)]
colnames(MG_Day30_Rep1_tpm)[2] <- "MG_Day30_Rep1"

# Claire's Day30 #1
MG_Day30_Rep2_genes <- read.table("MG_Day30_Rep2.merged.genes.results", sep="\t", header=TRUE)
MG_Day30_Rep2_genes[,1]<- unlist(strsplit(as.character(MG_Day30_Rep2_genes[, 1]), split="\\."))[(1:length(MG_Day30_Rep2_genes[, 1]))*2-1]
MG_Day30_Rep2_tpm <- MG_Day30_Rep2_genes[,c(1,6)]
colnames(MG_Day30_Rep2_tpm)[2] <- "MG_Day30_Rep2"

# Read in Glass data ####
# Ex_vivo #1
Microglia_ExVivo_Rep1_genes <- read.table("Microglia_ExVivo_Rep1.merged.genes.results", sep="\t", header=TRUE)
Microglia_ExVivo_Rep1_genes[,1]<- unlist(strsplit(as.character(Microglia_ExVivo_Rep1_genes[, 1]), split="\\."))[(1:length(Microglia_ExVivo_Rep1_genes[, 1]))*2-1]
Microglia_ExVivo_Rep1_tpm <- Microglia_ExVivo_Rep1_genes[,c(1,6)]
colnames(Microglia_ExVivo_Rep1_tpm)[2] <- "Microglia_ExVivo_Rep1"

# Ex_vivo #2
Microglia_ExVivo_Rep2_genes <- read.table("Microglia_ExVivo_Rep2.merged.genes.results", sep="\t", header=TRUE)
Microglia_ExVivo_Rep2_genes[,1]<- unlist(strsplit(as.character(Microglia_ExVivo_Rep2_genes[, 1]), split="\\."))[(1:length(Microglia_ExVivo_Rep2_genes[, 1]))*2-1]
Microglia_ExVivo_Rep2_tpm <- Microglia_ExVivo_Rep2_genes[,c(1,6)]
colnames(Microglia_ExVivo_Rep2_tpm)[2] <- "Microglia_ExVivo_Rep2"

# Ex_vivo #3
Microglia_ExVivo_Rep3_genes <- read.table("Microglia_ExVivo_Rep3.merged.genes.results", sep="\t", header=TRUE)
Microglia_ExVivo_Rep3_genes[,1]<- unlist(strsplit(as.character(Microglia_ExVivo_Rep3_genes[, 1]), split="\\."))[(1:length(Microglia_ExVivo_Rep3_genes[, 1]))*2-1]
Microglia_ExVivo_Rep3_tpm <- Microglia_ExVivo_Rep3_genes[,c(1,6)]
colnames(Microglia_ExVivo_Rep3_tpm)[2] <- "Microglia_ExVivo_Rep3"

# In_Vitro #1
Microglia_InVitro_Rep1_genes <- read.table("Microglia_InVitro_Rep1.merged.genes.results", sep="\t", header=TRUE)
Microglia_InVitro_Rep1_genes[,1]<- unlist(strsplit(as.character(Microglia_InVitro_Rep1_genes[, 1]), split="\\."))[(1:length(Microglia_InVitro_Rep1_genes[, 1]))*2-1]
Microglia_InVitro_Rep1_tpm <- Microglia_InVitro_Rep1_genes[,c(1,6)]
colnames(Microglia_InVitro_Rep1_tpm)[2] <- "Microglia_InVitro_Rep1"

# In_Vitro #2
Microglia_InVitro_Rep2_genes <- read.table("Microglia_InVitro_Rep2.merged.genes.results", sep="\t", header=TRUE)
Microglia_InVitro_Rep2_genes[,1]<- unlist(strsplit(as.character(Microglia_InVitro_Rep2_genes[, 1]), split="\\."))[(1:length(Microglia_InVitro_Rep2_genes[, 1]))*2-1]
Microglia_InVitro_Rep2_tpm <- Microglia_InVitro_Rep2_genes[,c(1,6)]
colnames(Microglia_InVitro_Rep2_tpm)[2] <- "Microglia_InVitro_Rep2"

# Read in iPSC ####


# Rep 1
iPSCs_Rep1 <- read.table("iPSCs_Rep1.merged.genes.results", sep="\t", header=TRUE)
iPSCs_Rep1[,1]<- unlist(strsplit(as.character(iPSCs_Rep1[, 1]), split="\\."))[(1:length(iPSCs_Rep1[, 1]))*2-1]
iPSCs_Rep1_tpm <- iPSCs_Rep1[,c(1,6)]
colnames(iPSCs_Rep1_tpm)[2] <- "iPSC_1"

# Rep 2
iPSCs_Rep2 <- read.table("iPSCs_Rep2.merged.genes.results", sep="\t", header=TRUE)
iPSCs_Rep2[,1]<- unlist(strsplit(as.character(iPSCs_Rep2[, 1]), split="\\."))[(1:length(iPSCs_Rep2[, 1]))*2-1]
iPSCs_Rep2_tpm <- iPSCs_Rep2[,c(1,6)]
colnames(iPSCs_Rep2_tpm)[2] <- "iPSC_2"

# Read in the Monocyte Data ####


# Add Monocyte derived Macrophage rep1
Mon_derived_Macro_1_genes <- read.table("Monocyte_Rep1.merged.genes.results", sep="\t", header=TRUE)
Mon_derived_Macro_1_genes[,1]<- unlist(strsplit(as.character(Mon_derived_Macro_1_genes[, 1]), split="\\."))[(1:length(Mon_derived_Macro_1_genes[, 1]))*2-1]
Mon_derived_Macro_1_tpm <- Mon_derived_Macro_1_genes[,c(1,6)]
colnames(Mon_derived_Macro_1_tpm)[2] <- "Mon_derived_Macro_1"

# Add Monocyte derived Macrophage rep1
Mon_derived_Macro_2_genes <- read.table("Monocyte_Rep2.merged.genes.results", sep="\t", header=TRUE)
Mon_derived_Macro_2_genes[,1]<- unlist(strsplit(as.character(Mon_derived_Macro_2_genes[, 1]), split="\\."))[(1:length(Mon_derived_Macro_2_genes[, 1]))*2-1]
Mon_derived_Macro_2_tpm <- Mon_derived_Macro_2_genes[,c(1,6)]
colnames(Mon_derived_Macro_2_tpm)[2] <- "Mon_derived_Macro_2"


# Add Hum_blood_derived_Macro rep1
Hum_blood_derived_Macro_1_genes <- read.table("Macrophage_Rep1.merged.genes.results", sep="\t", header=TRUE)
Hum_blood_derived_Macro_1_genes[,1]<- unlist(strsplit(as.character(Hum_blood_derived_Macro_1_genes[, 1]), split="\\."))[(1:length(Hum_blood_derived_Macro_1_genes[, 1]))*2-1]
Hum_blood_derived_Macro_1_tpm <- Hum_blood_derived_Macro_1_genes[,c(1,6)]
colnames(Hum_blood_derived_Macro_1_tpm)[2] <- "Hum_blood_derived_Macro_1"

# AddHum_blood_derived_Macrorep1
Hum_blood_derived_Macro_2_genes <- read.table("Macrophage_Rep2.merged.genes.results", sep="\t", header=TRUE)
Hum_blood_derived_Macro_2_genes[,1]<- unlist(strsplit(as.character(Hum_blood_derived_Macro_2_genes[, 1]), split="\\."))[(1:length(Hum_blood_derived_Macro_2_genes[, 1]))*2-1]
Hum_blood_derived_Macro_2_tpm <- Hum_blood_derived_Macro_2_genes[,c(1,6)]
colnames(Hum_blood_derived_Macro_2_tpm)[2] <- "Hum_blood_derived_Macro_2"

# Combine TPM files (Need to work on making each gene unique and look into genes which are missing)
library(dplyr)
TPM <- inner_join(MG_genes,MG_Day2_Rep1_tpm, by="gene_id")
TPM <- inner_join(TPM,MG_Day2_Rep2_tpm, by="gene_id")
TPM <- inner_join(TPM,MG_Day7_Rep1_tpm, by="gene_id")
TPM <- inner_join(TPM,MG_Day7_Rep2_tpm, by="gene_id")
TPM <- inner_join(TPM,MG_Day14_Rep1_tpm, by="gene_id")
TPM <- inner_join(TPM,MG_Day14_Rep2_tpm, by="gene_id")
TPM <- inner_join(TPM,MG_Day30_Rep1_tpm, by="gene_id")
TPM <- inner_join(TPM,MG_Day30_Rep2_tpm, by="gene_id")
TPM <- inner_join(TPM,Microglia_ExVivo_Rep1_tpm, by="gene_id")
TPM <- inner_join(TPM,Microglia_ExVivo_Rep2_tpm, by="gene_id")
TPM <- inner_join(TPM,Microglia_ExVivo_Rep3_tpm, by="gene_id")
TPM <- inner_join(TPM,Microglia_InVitro_Rep1_tpm, by="gene_id")
TPM <- inner_join(TPM,Microglia_InVitro_Rep2_tpm, by="gene_id")
TPM <- inner_join(TPM,iPSCs_Rep1_tpm, by="gene_id")
TPM <- inner_join(TPM,iPSCs_Rep2_tpm, by="gene_id")
TPM <- inner_join(TPM,Mon_derived_Macro_1_tpm, by="gene_id")
TPM <- inner_join(TPM,Mon_derived_Macro_2_tpm, by="gene_id")
TPM <- inner_join(TPM,Hum_blood_derived_Macro_1_tpm, by="gene_id")
TPM <- inner_join(TPM,Hum_blood_derived_Macro_2_tpm, by="gene_id")
TPM <- unique(TPM)

# Make Final table
TPM_Final <- data.frame(TPM[,1])
row.names(TPM_Final) <- TPM[,1]
TPM_Final<-TPM_Final[,-1]

# Take the mean of Claire's Microglia Day 2
TPM_Final$MG_Day2 <- rowMeans(TPM[,c('MG_Day2_Rep1','MG_Day2_Rep2')])

# Take the mean of Claire's Microglia Day 7
TPM_Final$MG_Day7 <- rowMeans(TPM[,c('MG_Day7_Rep1','MG_Day7_Rep2')])
                                  
# Take the mean of Claire's Microglia Day 14
TPM_Final$MG_Day14 <- rowMeans(TPM[,c('MG_Day14_Rep1','MG_Day14_Rep2')]) 

# Take the mean of Claire's Microglia Day 30
TPM_Final$MG_Day30 <- rowMeans(TPM[,c('MG_Day30_Rep1','MG_Day30_Rep2')])

# Take the mean of Glass ExVivo
TPM_Final$Glass_ExVivo <- rowMeans(TPM[,c("Microglia_ExVivo_Rep1","Microglia_ExVivo_Rep2","Microglia_ExVivo_Rep3")])

# Take the mean of Glass InVivo
TPM_Final$Glass_InVitro <- rowMeans(TPM[,c("Microglia_InVitro_Rep1","Microglia_InVitro_Rep2")])

# Take the mean of iPSC
TPM_Final$iPSC <- rowMeans(TPM[,c('iPSC_1','iPSC_2')])

# Take the mean of Mon_derived_Macro
TPM_Final$Mono_derived_Macro <- rowMeans(TPM[,c("Mon_derived_Macro_1","Mon_derived_Macro_2")])

# Take the mean of Hum_blood_derived_Macro
TPM_Final$Hum_blood_derived_Macro <- rowMeans(TPM[,c("Hum_blood_derived_Macro_1","Hum_blood_derived_Macro_2")])

# Change the names
colnames(TPM_Final) <- c('GFiMG Day2', 'GFiMG Day7', 'GFiMG Day14', 'GFiMG Day30', 'Ex vivo microglia (Gosselin et al. 2017)', 'In vitro Day7 microglia (Gosselin et al. 2017)','WTC iPSC','monocyte-derived macrophage (Heinz et al, 2018)','monocyte-derived macrophage (Carlin et al, 2018)')

# perform Log transformation
TPM_Final_Log2 <- log2(TPM_Final+1)
TPM_Final_Log10 <- log10(TPM_Final+1)

# Make Heat Map
library(pheatmap)

png('RNA-seq_log2_pheatmap.png',res=200,units="in",width=8,height=12)
pheatmap(TPM_Final_Log2,main = "849 Highly Expressed Genes (Log2 TPM)",show_rownames = F, cluster_rows = T, cluster_cols = T, clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean")
dev.off()

png('RNA-seq_log10_pheatmap.png',res=200,units="in",width=8,height=12)
pheatmap(TPM_Final_Log10, main = "849 Highly Expressed Genes (Log10 TPM)" , show_rownames = F, labels_row = NA, cluster_rows = T, cluster_cols = T, clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean")
dev.off()

