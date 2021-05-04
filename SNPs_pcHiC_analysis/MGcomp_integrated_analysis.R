library(rtracklayer)
library(GenomicRanges)
library(genomation) 
load('/shen/shenlabstore3/neng/20200309_MGEN/SNPs_pcHiC_analysis/data/MGcomp/MGcomp.atac.seq.peaks.ann.Rdata')
load('/shen/shenlabstore3/neng/20200309_MGEN/SNPs_pcHiC_analysis/data/MGcomp/MGcomp.atac.seq.peaks.Rdata') 


celltypes<-c('Microglia','Microglia_InVitro', 'Microglia_ExVivo')
initials<-c('MG','MG_InVitro','MG_ExVivo')



getdistalATAC<- function(celltype, ann, res, initial) {
	distalID<-which(ann$'peak_type' == "distal")
	distalATAC<-res[distalID,]
	bedFile<-makeGRangesFromDataFrame(distalATAC)
	outbed<-resize(bedFile, width=500, fix="center")
	outfile<-paste('final_MGcomp/',initial,'.distalATAC.500bp.bed',sep="")
	export.bed(outbed, con=outfile)
	outfileproc<-paste('final_MGcomp/',initial,'.distalATAC.500bp.srt.merged.bed',sep="")
	processfiles(outfile,outfileproc)
}

processfiles <- function(file, ofile) {
	command<-paste('sort -k1,1 -k2,2n',file ,'| bedtools merge -i - >', ofile,sep=' ')
	system(command)
}

for (i in 1:length(celltypes)){
	celltype<-celltypes[i]
	#print(celltype)
	initial<-initials[i]
	atac_ann<-atac.seq.peaks.ann[[celltype]]
	atac_res<-atac.seq.peaks[[celltype]]
	getdistalATAC(celltype, atac_ann, atac_res, initial)
}