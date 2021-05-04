library(DESeq2)
library(stringr)
library(ggplot2) 
library(RColorBrewer)
library(pheatmap)
library(edgeR)


runDESeq2 <- function(count.table, samples, condition) {

	dds <- DESeqDataSetFromMatrix(countData = count.table,colData = samples,design = ~ condition)

	keep<-rowSums(counts(dds)) > 1
	dds<- dds[keep,]
	dds <- DESeq(dds)
	res<-results(dds)
	resLFC <- lfcShrink(dds, coef=2, type="apeglm")

	#save(res, file='/shen/shenlabstore3/neng/20200617_NaiveForm/Interaction_Analysis/featureCounts/deseq2.results.Rdata')

	vsd <- vst(dds, blind=FALSE)
	rld <- rlog(dds, blind=FALSE)
	sampleDists <- dist(t(assay(vsd)))
	DESeq2::plotMA(resLFC, ylim=c(-4,4))


	#pdf('RNA_heatmap.pdf',height=7,width=7)
	#sampleDistMatrix <- as.matrix(sampleDists)
	#rownames(sampleDistMatrix) <- samples.names
	#colnames(sampleDistMatrix) <- NULL
	#colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
	#pheatmap(sampleDistMatrix,
	#         clustering_distance_rows=sampleDists,
	#         clustering_distance_cols=sampleDists,
	#         col=colors)

	#dev.off()


	#pdf('RNA_pca.pdf',height=7,width=7)
	#plotPCA(rld)
	#dev.off()
	return(res)


}




summarizeExpressionResults <- function(expression.results, gene.lengths) {
	# Create results matrix where rows are genes and columns are cell types.
	cell.types <- names(expression.results)
	num.genes <- length(expression.results[[cell.types[1]]][, 1])
	expression.summary <- data.frame(matrix(NA, num.genes, length(cell.types)), stringsAsFactors=F)
	colnames(expression.summary) <- cell.types
	rownames(expression.summary) <- rownames(expression.results[[cell.types[1]]])

	# Count number of replicates across all cell types.
	num.replicates <- 0
	for (cell.type in cell.types)
	  num.replicates <- num.replicates + length(names(expression.results[[cell.type]]))

	# Compile matrices containing gene lengths and expected counts for each replicate across all cell types.
	all.results <- data.frame(matrix(NA, num.genes, num.replicates))
	rownames(all.results) <- rownames(expression.results[[cell.types[1]]])

	groups <- c()
	col <- 1
	for (cell.type in cell.types) {
	  
	  all.results[, col:(col + length(names(expression.results[[cell.type]])) - 1)] <- expression.results[[cell.type]]
	  colnames(all.results)[col:(col + length(names(expression.results[[cell.type]])) - 1)] <- names(expression.results[[cell.type]])
	  groups <- c(groups, rep(cell.type, length(expression.results[[cell.type]])))
	  col <- col + length(names(expression.results[[cell.type]]))
	  
	}

	# Calculate TMM normalization factors for each replicate using edgeR.
	y <- DGEList(counts=all.results, group=groups)
	y <- calcNormFactors(y, method="TMM")

	# Calculate the TMM-normalized RPKM for each replicate (accounting for its reported gene lengths) using edgeR.
	all.rpkms <- data.frame(matrix(NA, num.genes, num.replicates))
	for (i in 1:length(all.results)) {
	  
	  RPKM <- rpkm(y, gene.length=gene.lengths[,1])
	  all.rpkms[, i] <- RPKM[, i]
	  
	}

	# Calculate the mean TMM-normalized RPKM across all replicates.
	for (i in 1:length(cell.types)) {
	  expression.summary[, i] <- rowMeans(as.matrix(all.rpkms[, which(groups==cell.types[i])]))
	}
	return(expression.summary)
}


counts<-read.csv("RNA_merged_counts.txt", sep="", head=T, skip=1, row.names = "Geneid")
geneids<-rownames(counts)
geneids.fixed<-unlist(lapply(geneids,function(gene) gsub('\\..*','',gene)))
rownames(counts)<-geneids.fixed
samples.names<-str_split_fixed(str_split_fixed(colnames(counts)[6:11],'output.',2)[,2],'.Aligned.',2)[,1]
colnames(counts)[6:11]<-samples.names


#### DESeq Analysis ####
condition<-str_split_fixed(samples.names,'_',2)[,1]
samples<-cbind(colnames(counts)[6:11],condition)
rownames(samples)<-samples[,1]
samples<-as.data.frame(samples[,-1])
colnames(samples)<-"condition"
count.table<-counts[,6:11]
genes<-rownames(counts)


deseq2.results<-runDESeq2(count.table, samples, condition)

expression.results<-list()
expression.results[['naive']]<-count.table[1:3]
expression.results[['form']]<-count.table[4:6]
gene.lengths<-counts[5]
normalized.expression<-summarizeExpressionResults(expression.results,gene.lengths)
normalized.expression<-normalized.expression[order(row.names(normalized.expression)),]

write.table(normalized.expression, file="/shen/shenlabstore3/neng/20200617_NaiveForm/Interaction_Analysis/featureCounts/
	normalized_expression.txt", sep="\t", row.names=T, col.names=T, quote=F)
write.table(deseq2.results, file="/shen/shenlabstore3/neng/20200617_NaiveForm/Interaction_Analysis/featureCounts/deseq2_results.txt", sep="\t", row.names=T, col.names=T, quote=F)

save(normalized.expression,file="/shen/shenlabstore3/neng/20200617_NaiveForm/Interaction_Analysis/featureCounts/normalized.expression.Rdata")

save(deseq2.results,file="/shen/shenlabstore3/neng/20200617_NaiveForm/Interaction_Analysis/featureCounts/deseq2.results.Rdata")







