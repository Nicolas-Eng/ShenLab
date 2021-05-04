
library(rtracklayer)
library(GenomicRanges)
library(genomation) 
load('/shen/shenlabstore3/neng/20200309_MGEN/SNPs_pcHiC_analysis/data/MGEN/atac.seq.peaks.ann.Rdata')
load('/shen/shenlabstore3/neng/20200309_MGEN/SNPs_pcHiC_analysis/data/MGEN/atac.seq.peaks.Rdata')

load('/shen/shenlabstore3/neng/20200309_MGEN/SNPs_pcHiC_analysis/data/MGEN/interactions.sig.ann.Rdata')
load('/shen/shenlabstore3/neng/20200309_MGEN/SNPs_pcHiC_analysis/data/MGEN/interactions.sig.res.Rdata')


celltypes<-c('ExcitatoryNeurons','Microglia')
initials<-c('EN','MG')



getInteractionsxdistalgeneATAC <- function(celltype, ann, res, initial) {
	#LHS distal ATAC-seq peak
	ATAC_IDS_LHS<-which((ann$'distal_ATAC-seq_lhs' > 0) & (ann$'promoter_lhs' == 0) & (ann$'promoter_ATAC-seq_lhs' == 0)) 
	ATAC_IDS_RHS<-which((ann$'distal_ATAC-seq_rhs' > 0) & (ann$'promoter_rhs' == 0) & (ann$'promoter_ATAC-seq_rhs' == 0))
 	#RHS distal ATAC-seq peak
	Interactions_x_ATAC_IDS_LHS<-res[ATAC_IDS_LHS,][,1:3]
	Interactions_x_ATAC_IDS_RHS<-res[ATAC_IDS_RHS,][,5:7]
	file_LHS<-paste('tmp/', initial,'.Inter_x_ATACDistal.LHS.bed',sep="")
	file_RHS<-paste('tmp/', initial,'.Inter_x_ATACDistal.RHS.bed',sep="")
	ofile<-paste('final/',initial,'.Inter_x_ATACDistal.merge.sort.bed',sep="")
	writefiles(Interactions_x_ATAC_IDS_LHS,file_LHS)
	writefiles(Interactions_x_ATAC_IDS_RHS,file_RHS)
	processfiles(file_LHS,file_RHS,ofile)
}

getInteractionsxdistalgeneATACpreq <- function(celltype, ann, res, initial) {
	#RHS gene connection with LHS distal ATAC-seq peak
	ATAC_IDS_LHS<-which(ann$'distal_ATAC-seq_lhs' > 0 & (ann$'promoter_rhs' > 0 | ann$'promoter_ATAC-seq_rhs' > 0) & (ann$'promoter_lhs' == 0) & (ann$'promoter_ATAC-seq_lhs' == 0)) 
	#LHS gene connection with RHS distal ATAC-seq peak
	ATAC_IDS_RHS<-which(ann$'distal_ATAC-seq_rhs' > 0 & (ann$'promoter_lhs' > 0 | ann$'promoter_ATAC-seq_lhs' > 0) & (ann$'promoter_rhs' == 0) & (ann$'promoter_ATAC-seq_rhs' == 0)) 
	Interactions_x_ATAC_IDS_LHS<-res[ATAC_IDS_LHS,][,1:3]
	Interactions_x_ATAC_IDS_RHS<-res[ATAC_IDS_RHS,][,5:7]
	file_LHS<-paste('tmp/',initial,'.Inter_x_ATACDistal.LHS.preq.bed',sep="")
	file_RHS<-paste('tmp/',initial,'.Inter_x_ATACDistal.RHS.preq.bed',sep="")
	ofile<-paste('final/',initial,'.Inter_x_ATACDistal.merge.sort.preq.bed',sep="")
	writefiles(Interactions_x_ATAC_IDS_LHS,file_LHS)
	writefiles(Interactions_x_ATAC_IDS_RHS,file_RHS)
	processfiles(file_LHS,file_RHS,ofile)
}


getdistalInteractions <- function(celltype, ann, res, initial) {
	#LHS distal interaction connected to RHS gene
	InterDistal_IDS_LHS<- which((ann$'promoter_rhs' > 0 | ann$'promoter_ATAC-seq_rhs' > 0) & (ann$'promoter_lhs' == 0) & (ann$'promoter_ATAC-seq_lhs' == 0))  
	#RHS distal interaction connected to LHS gene
	InterDistal_IDS_RHS<- which((ann$'promoter_lhs' > 0 | ann$'promoter_ATAC-seq_lhs' > 0) & (ann$'promoter_rhs' == 0) & (ann$'promoter_ATAC-seq_rhs' == 0))
	InterDistal_LHS<-res[InterDistal_IDS_LHS,][,1:3]
	InterDistal_RHS<-res[InterDistal_IDS_RHS,][,5:7]
	file_LHS<-paste('tmp/',initial,'.InterDistal.LHS.preq.bed',sep="")
	file_RHS<-paste('tmp/',initial,'.InterDistal.RHS.preq.bed',sep="")
	ofile<-paste('final/',initial,'.InterDistal.merge.sort.preq.bed',sep="")
	writefiles(InterDistal_LHS,file_LHS)
	writefiles(InterDistal_RHS,file_RHS)
	processfiles(file_LHS,file_RHS,ofile)
}

writefiles <- function(var, file) {
	write.table(var,file, col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
}

processfiles <- function(file_LHS, file_RHS, ofile) {
	command<-paste('cat', file_LHS, file_RHS,'| sort -k1,1 -k2,2n - | bedtools merge -i - >', ofile,sep=' ')
	system(command)
}


getdistalATACexp<- function(celltype, ann, res, initial) {
	distalID<-which(ann$'peak_type' == "distal")
	distalATAC<-res[distalID,]
	bedFile<-makeGRangesFromDataFrame(distalATAC)
	outbed<-resize(bedFile, width=1000, fix="center")
	export.bed(outbed, con=paste('final/',initial,'.distalATAC.1kb.bed',sep=""))
}

getdistalATAC<- function(celltype, ann, res, initial) {
	distalID<-which(ann$'peak_type' == "distal")
	distalATAC<-res[distalID,]
	bedFile<-makeGRangesFromDataFrame(distalATAC)
	ifile<-paste('final/',initial,'.distalATAC.bed',sep="")
	export.bed(bedFile, con=ifile)
	distalinteract<-paste('final/',initial,'.InterDistal.merge.sort.preq.bed',sep="")
	ofile<-paste('final/',initial,'.distalATACxInterDistal.merge.sort.preq.bed',sep="")
	command<-paste('bedtools intersect -a', ifile, '-b', distalinteract, '>', ofile, sep=" ")
	system(command)
	obedFile<-readBed(ofile)
	outbed<-resize(obedFile,width=1000, fix="center")
	export.bed(outbed, con=paste('final/',initial,'.distalATACxInterDistal.merge.sort.preq.1kb.bed',sep=""))
}




for (i in 1:length(celltypes)){
	celltype<-celltypes[i]
	#print(celltype)
	initial<-initials[i]
	interact_ann<-interactions.sig.ann[[celltype]]
	interact_res<-interactions.sig.res[[celltype]]
	#getInteractionsxdistalgeneATAC(celltype, interact_ann, interact_res, initial)
	#getInteractionsxdistalgeneATACpreq(celltype, interact_ann, interact_res, initial)
	#getdistalInteractions(celltype, interact_ann, interact_res, initial)
	atac_ann<-atac.seq.peaks.ann[[celltype]]
	atac_res<-atac.seq.peaks[[celltype]]
	getdistalATAC(celltype, atac_ann, atac_res, initial)
	#getdistalATACexp(celltype, atac_ann, atac_res, initial)
}








