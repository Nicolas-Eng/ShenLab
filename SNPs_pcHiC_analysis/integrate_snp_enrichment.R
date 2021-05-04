library(rtracklayer)
library(GenomicRanges)
library(genomation) 

load('/shen/shenlabstore3/neng/20200309_MGEN/SNPs_pcHiC_analysis/data/MGEN/interactions.sig.ann.Rdata')
load('/shen/shenlabstore3/neng/20200309_MGEN/SNPs_pcHiC_analysis/data/MGEN/interactions.sig.res.Rdata')


celltypes<-c('ExcitatoryNeurons')
initials<-c('EN')

getdistalInteractions <- function(celltype, ann, res, initial) {
	#LHS distal interaction connected to RHS gene
	InterDistal_IDS_LHS<- which((ann$'promoter_rhs' > 0 | ann$'promoter_ATAC-seq_rhs' > 0) & (ann$'promoter_lhs' == 0) & (ann$'promoter_ATAC-seq_lhs' == 0))  
	#RHS distal interaction connected to LHS gene
	InterDistal_IDS_RHS<- which((ann$'promoter_lhs' > 0 | ann$'promoter_ATAC-seq_lhs' > 0) & (ann$'promoter_rhs' == 0) & (ann$'promoter_ATAC-seq_rhs' == 0))

	x<-res[InterDistal_IDS_LHS,]
	y<-res[InterDistal_IDS_RHS,]
	out<-rbind(x,y)
	print(out)
	save(out,'DitalInteracting')
#	file_LHS<-paste('tmp_SNPs_enrich/',initial,'.InterDistal.LHS.full.preq.bed',sep="")
#	file_RHS<-paste('tmp_SNPs_enrich/',initial,'.InterDistal.RHS.full.preq.bed',sep="")
#	ofile<-paste('final_SNPs_enrich/',initial,'.InterDistal.merged.full.preq.bed',sep="")
#	writefiles(InterDistal_LHS_corrected,file_LHS)
#	writefiles(InterDistal_RHS,file_RHS)
#	processfiles(file_LHS,file_RHS,ofile)
}


writefiles <- function(var, file) {
	write.table(var,file, col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
}

processfiles <- function(file_LHS, file_RHS, ofile) {
	command<-paste('cat', file_LHS, file_RHS,'>', ofile,sep=' ')
	system(command)
}


for (i in 1:length(celltypes)){
	celltype<-celltypes[i]
	#print(celltype)
	initial<-initials[i]
	interact_ann<-interactions.sig.ann[[celltype]]
	interact_res<-interactions.sig.res[[celltype]]
	getdistalInteractions(celltype, interact_ann, interact_res, initial)
}


