library(bedr)
library(data.table)
library(plyr)

load('files/XOR.interactions.annotated.Rdata')
#Correct the annotation so that it is based on gene names rather than interaction in 5kb bins

cell.types<-c('mESCs','epiLCs')

correctInteractions<-function(cell.types, XOR.interactions.annotated) {
	XOR.interactions.gene.corrected <-list()
	for (cell.type in cell.types) {
		df.corrected <- data.frame(anchorcoord=character(),
	                distalcoord=character(), 
	                fdr=character(),
	                gene=character(), 
	                stringsAsFactors=FALSE)
		for (i in 1:length(XOR.interactions.annotated[[cell.type]][,1])) {
			print(i)
			gene.split<-unlist(strsplit(XOR.interactions.annotated[[cell.type]][['gene']][i],','))
			anchorcoord<-XOR.interactions.annotated[[cell.type]][i,]['anchorcoords'][[1]]
			distalcoord<-XOR.interactions.annotated[[cell.type]][i,]['distalcoords'][[1]]
			fdr<-XOR.interactions.annotated[[cell.type]][i,]['fdr'][[1]]
			for (gene in gene.split) {
				gene<-gsub('\\..*','',gene)
				newrow<-data.frame(anchorcoord,distalcoord,fdr,gene,stringsAsFactors=FALSE)
				df.corrected<-rbind(df.corrected,newrow)
				print(newrow)
			}
		}
		XOR.interactions.gene.corrected[[cell.type]]<-df.corrected
	}
	return(XOR.interactions.gene.corrected)
}


collapseAnnotations <- function(df,varColBy, annotate.list ) {
	abedr.dt.tmp<-data.table(df)
	abedr.dt.out<-data.table()
	for (ann in annotate.list){
		abedr.dt.collapsed<-abedr.dt.tmp[,paste0(get(ann),collapse=","),by=c(as.character(varColBy))]
		abedr.dt.collapsed$V1<-unlist(lapply(abedr.dt.collapsed$V1, function(x) gsub(" ","",toString(unique(unlist(strsplit(x,",")))))))
		colnames(abedr.dt.collapsed)<-c(varColBy,ann)
		#print(head(abedr.dt.collapsed))
		if (empty(abedr.dt.out)) {
			abedr.dt.out <- abedr.dt.collapsed
		} else {
			abedr.dt.out <- merge(abedr.dt.out, abedr.dt.collapsed,by=varColBy)
		}
	}
	print(head(abedr.dt.out),2)
	abedr.df.out<-setDF(abedr.dt.out)
	abedr.df.out[[1]]<-as.character(abedr.df.out[[1]])
	return(abedr.df.out)
}

#XOR.interactions.gene.corrected<-correctInteractions(cell.types, XOR.interactions.annotated)
#save(XOR.interactions.gene.corrected,file='files/XOR.interactions.gene.corrected.Rdata')

load('files/XOR.interactions.gene.corrected.Rdata')

RSEM.table<-read.table('RSEM/rpkm.data.individual.txt',stringsAsFactors=FALSE,header=TRUE)


#Get rid of genes with < 1 RPKM for both cell-types
RSEM.df.processed<-RSEM.table[!(RSEM.table[['mESCs']]<=1 & RSEM.table[['epiLCs']]<=1),]

category <- vector(mode = "list", length = length(RSEM.df.processed[,1]))
for (i in 1:length(RSEM.df.processed[,1])) {
	if (RSEM.df.processed[['mESCs']][i] < RSEM.df.processed[['epiLCs']][i]) {
		category[[i]]<-'increase'
	} else if (RSEM.df.processed[['epiLCs']][i] < RSEM.df.processed[['mESCs']][i]) {
		category[[i]]<-'decrease'
	} else {
		category[[i]]<-'equals'
	}
}

RSEM.df.processed.category<-cbind(RSEM.df.processed,data.frame(unlist(category)))



sumFDRinteractions <- function(cell.types, XOR.interactions.gene.corrected) {
	XOR.interactions.gene.fdr.corrected <-list()
	for (cell.type in cell.types) {
		XOR.cell.type <- XOR.interactions.gene.corrected[[cell.type]]
		XOR.gene.collapsed.fdr<-collapseAnnotations(XOR.cell.type, "gene", "fdr")
		XOR.gene.collapsed.anchor<-collapseAnnotations(XOR.cell.type, "gene", "anchorcoord")
		XOR.gene.collapsed.distal<-collapseAnnotations(XOR.cell.type, "gene", "distalcoord")

		XOR.gene.collapsed.comb<-join(XOR.gene.collapsed.anchor,XOR.gene.collapsed.distal,'gene')
		XOR.gene.collapsed.comb<-join(XOR.gene.collapsed.comb,XOR.gene.collapsed.fdr)
		XOR.gene.collapsed.comb<-XOR.gene.collapsed.comb[c("anchorcoord","distalcoord","fdr","gene")]

		XOR.gene.collapsed.comb[['fdr']]<-unlist(lapply(XOR.gene.collapsed.comb[['fdr']],
			function(x) sum(-log(as.numeric(unlist((strsplit(unlist(x),","))))))))
		XOR.interactions.gene.fdr.corrected[[cell.type]]<-XOR.gene.collapsed.comb
	}
	return(XOR.interactions.gene.fdr.corrected)
}


XOR.interactions.gene.fdr.corrected<-sumFDRinteractions(cell.types,XOR.interactions.gene.corrected)
XOR.genes<-unique(c(XOR.interactions.gene.fdr.corrected[['mESCs']][['gene']],XOR.interactions.gene.fdr.corrected[['epiLCs']][['gene']]))
XOR.genes.trends<-na.omit(RSEM.df.processed.category[XOR.genes,])
colnames(XOR.genes.trends)<-c('RPKM.mESCs','RPKM.epiLCs','RPKM.trend')
genes<-rownames(XOR.genes.trends)



load('RSEM/DESeq.Rdata')
DEseq.na.omit<-na.omit(resOrdered)
DEseq.genes<- rownames(DEseq.na.omit)
DE.genes<-DEseq.genes[DEseq.na.omit$padj < 0.000001]


mESCs.fdr <- vector(mode = "list", length = length(genes))
epiLCs.fdr <- vector(mode = "list", length = length(genes))
interaction.trend <- vector(mode = "list", length = length(genes))
DE.trend <- vector(mode="list", length=length(genes))
i=1
for (gene in genes) {
	if (gene %in% XOR.interactions.gene.fdr.corrected[['mESCs']][['gene']]) {
		mESCs.index<-which(XOR.interactions.gene.fdr.corrected[['mESCs']][['gene']]==gene)
		mESCs.fdr[[i]]<-XOR.interactions.gene.fdr.corrected[['mESCs']][mESCs.index,][['fdr']]
	}
	else {
		mESCs.fdr[[i]]<-0
	}
	mESCs.fdr
	if (gene %in% XOR.interactions.gene.fdr.corrected[['epiLCs']][['gene']]) {
		epiLCs.index<-which(XOR.interactions.gene.fdr.corrected[['epiLCs']][['gene']]==gene)
		epiLCs.fdr[[i]]<-XOR.interactions.gene.fdr.corrected[['epiLCs']][epiLCs.index,][['fdr']]
	}
	else {
		epiLCs.fdr[[i]]<-0
	}
	if (epiLCs.fdr[[i]] < mESCs.fdr[[i]]) {
		interaction.trend[[i]]<-'decrease'
	} else if (epiLCs.fdr[[i]] > mESCs.fdr[[i]]) {
		interaction.trend[[i]]<-'increase'
	} else {
		interaction.trend[[i]]<-'equals'
	}
	if (gene %in% DE.genes) {
		DE.trend[[i]] <- 'yes'
	} else {
		DE.trend[[i]] <- 'no'
	}
	i=i+1
}

epiLCs.fdr<-unlist(epiLCs.fdr)
mESCs.fdr<-unlist(mESCs.fdr)
DE.trend<-unlist(DE.trend)
interaction.trend<-unlist(interaction.trend)
interaction.trend.df<-data.frame(cbind(genes,mESCs.fdr,epiLCs.fdr,interaction.trend,DE.trend),stringsAsFactors=FALSE)
colnames(interaction.trend.df)<-c('gene','mESCs.fdr','epiLCs.fdr','interaction.trend','DEseq')


genes.trend.df<-data.frame((setDT(XOR.genes.trends,keep.rownames=TRUE)[]),stringsAsFactors=FALSE)
colnames(genes.trend.df)<-c('gene','mESCs.RPKM','epiLCs.RPKM','gene.trend')

trend.merged.df<-join(interaction.trend.df, genes.trend.df,'gene')
#trend.merged.df.increase<-trend.merged.df[trend.merged.df[['interaction.trend']]=='increase',]

#trend.merged.df.decrease<-trend.merged.df[trend.merged.df[['interaction.trend']]=='decrease',]

#diff.MAPS.increase<-as.numeric(trend.merged.df.increase[['epiLCs.fdr']])-as.numeric(trend.merged.df.increase[['mESCs.fdr']])
#diff.RPKM.increase<-log2(trend.merged.df.increase[['epiLCs.RPKM']]+1) - log2(trend.merged.df.increase[['mESCs.RPKM']]+1)

#diff.MAPS.decrease<-as.numeric(trend.merged.df.decrease[['epiLCs.fdr']])-as.numeric(trend.merged.df.decrease[['mESCs.fdr']])
#diff.RPKM.decrease<-log2(trend.merged.df.decrease[['epiLCs.RPKM']]+1) - log2(trend.merged.df.decrease[['mESCs.RPKM']]+1)

#cor.value.increase <- round(cor(diff.MAPS.increase, diff.RPKM.increase),4)
#cor.p.increase <- formatC(cor.test(diff.MAPS.increase, diff.RPKM.increase)$p.value, format = "e", digits = 2)

#cor.value.decrease <- round(cor(diff.MAPS.decrease, diff.RPKM.decrease),4)
#cor.p.decrease <- formatC(cor.test(diff.MAPS.decrease, diff.RPKM.decrease)$p.value, format = "e", digits = 2)

#plot( diff.MAPS.decrease, diff.RPKM.decrease, xlim=c(-60,60), ylim=c(-3,3), cex = 0.2, pch=19, col='gray')
#points( diff.MAPS.increase, diff.RPKM.increase, xlim=c(-60,60), ylim=c(-3,3), cex = 0.2, pch=19, col='blue')
#smoothScatter( diff.MAPS.increase, diff.RPKM.increase, xlim=c(-60,60), ylim=c(-7,7), cex = 0.2, pch=19, col='red')
#smoothScatter( diff.MAPS.decrease, diff.RPKM.decrease, xlim=c(-60,60), ylim=c(-7,7), cex = 0.1, pch=19, col='blue')
#xlab=paste('Difference in MAPS interactions: epiLCs - mESCs')
#ylab=paste('Difference in gene expression Log2(RPKM+1): epiLCs - mESCs')


#diff.MAPS<-as.numeric(trend.merged.df[['epiLCs.fdr']])-as.numeric(trend.merged.df[['mESCs.fdr']])
#diff.RPKM<-log2(trend.merged.df[['epiLCs.RPKM']]+1) - log2(trend.merged.df[['mESCs.RPKM']]+1)
#cor.value <- round(cor(diff.MAPS,diff.RPKM),4)


trend.increase.decrease.df<-trend.merged.df[trend.merged.df[['interaction.trend']]=='increase' & 
	trend.merged.df[['gene.trend']]=='decrease',]

trend.increase.increase.df<-trend.merged.df[trend.merged.df[['interaction.trend']]=='increase' & 
	trend.merged.df[['gene.trend']]=='increase',]

trend.decrease.decrease.df<-trend.merged.df[trend.merged.df[['interaction.trend']]=='decrease' &
	trend.merged.df[['gene.trend']]=='decrease',]

trend.decrease.increase.df<-trend.merged.df[trend.merged.df[['interaction.trend']]=='decrease' & 
	trend.merged.df[['gene.trend']]=='increase',]

save(trend.merged.df,file='downstream_new_0.000001_test/trend.merged.df.Rdata')


