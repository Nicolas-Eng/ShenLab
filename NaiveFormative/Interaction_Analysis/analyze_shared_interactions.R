library(bedr)
library(data.table)
library(plyr)

options(scipen=999)
#Format bed as chr:start-end from 3 columns of chr\tstart\tend
getBedrFormatted <- function(bed) {
	bedr <- paste0(bed[,1],':',bed[,2],'-',bed[,3])
	return(bedr)
}

#Format the RHS coordinates with bedr formatting chr:start-end
regroupAnnotationFile <- function(bed_inter) {
	tmpbed <- cbind(bed_inter[,2], bed_inter[,3], bed_inter[,4])
	#print(tmpbed)
	tmpbedr<-getBedrFormatted(tmpbed)
	bedr_inter <- cbind(bed_inter[[1]],tmpbedr)

	return(bedr_inter)
}

load('files/XOR.shared.switch.gene.interactions.Rdata')

cell.types<-c('mESCs','epiLCs')
XOR.switch<-c()
#intersect.atac.x.switch<-c()
for (cell.type in cell.types) {
	XOR.switch.gene<-as.character(XOR.shared.switch.gene.interactions[['gene']])
	XOR.switch.interactions<-as.character(XOR.shared.switch.gene.interactions[[paste0(cell.type,'.only')]])
	XOR.switch.interactions.coords<-unlist(lapply(XOR.switch.interactions,function(x) unlist((strsplit(x,",")))))
	XOR.switch.interactions.distal.only<-unique(bedr.sort.region(unlist(lapply(XOR.switch.interactions.coords,function(x)unlist(strsplit(x,'_'))[2]))))
#	bed.atac<-read.table(paste0('peaks/atac/',cell.type,'.optimal.atac.bed'),sep="\t", header=F, stringsAsFactors=F)
#	bedr.atac<-bedr.sort.region(getBedrFormatted(bed.atac))
#	intersect.atac.XOR.switch.interactions <- bedr(input = list(a=bedr.atac, b=XOR.switch.interactions.distal.only), method = "intersect", params = "-wo")
#	intersect.atac.x.switch[[cell.type]]<-intersect.atac.XOR.switch.interactions
#	XOR.switch.df<-convert2bed(XOR.switch.interactions.distal.only)
#	write.table(XOR.switch.df,paste0('downstream/',cell.type,'.switch.distal.bed'),quote=F,col.names=F,sep="\t",row.names=F)
#	atac.x.switch.df<-convert2bed(intersect.atac.x.switch[[cell.type]][['index']])
#	write.table(atac.x.switch.df,paste0('downstream/',cell.type,'.atac.x.switch.bed'),quote=F,col.names=F,sep="\t",row.names=F)

#	intersect.atac.x.switch[[cell.type]]<-regroupAnnotationFile(intersect.atac.x.switch[[cell.type]])
}


genes.switch<-unlist(lapply(as.character(XOR.shared.switch.gene.interactions[['gene']]),function(x) strsplit(x,",")))
genes.switch<-gsub("\\..*","",genes.switch)

rpkm.table<-read.table('RSEM/rpkm.data.individual.txt',sep="\t",header=T,stringsAsFactors=F)
rpkm.table[['logFC']]<-log2((rpkm.table[['epiLCs']]+1)/(rpkm.table[['mESCs']]+1))


rpkm.genes.switch<-rpkm.table[genes.switch,]
plot(log2(rpkm.genes.switch$mESCs+1),log2(rpkm.genes.switch$epiLCs+1),main='Enhancer Switch Genes log(RPKM)',xlab='mESCs',ylab='epiLCs')
cor(log2(rpkm.genes.switch$mESCs+1),log2(rpkm.genes.switch$epiLCs+1))
boxplot(log2(rpkm.genes.switch$mESCs+1),log2(rpkm.genes.switch$epiLCs+1))

load('files/XOR.mESCs.distinct.gene.interactions.Rdata')
genes.mESCs<-unlist(lapply(as.character(XOR.mESCs.distinct.gene.interactions[['gene']]),function(x) strsplit(x,",")))
genes.mESCs<-gsub("\\..*","",genes.mESCs)
genes.mESCs.corrected<-setdiff(genes.mESCs,genes.switch)
rpkm.genes.mESCs<-rpkm.table[genes.mESCs.corrected,]
plot(log2(rpkm.genes.mESCs$mESCs+1),log2(rpkm.genes.mESCs$epiLCs+1))
cor(log2(rpkm.genes.mESCs$mESCs+1),log2(rpkm.genes.mESCs$epiLCs+1))

mESCs.distinct.coords<-unlist(lapply(as.character(XOR.mESCs.distinct.gene.interactions[['XOR.mESCs.coords']]),function(x) strsplit(x,",")))
mESCs.distinct.distal.only<-unique(bedr.sort.region(unlist(lapply(mESCs.distinct.coords,function(x) unlist(strsplit(x,"_"))[2]))))
mESCs.distinct.df<-convert2bed(mESCs.distinct.distal.only)
write.table(mESCs.distinct.df,paste0('downstream/mESCs.distinct.distal.bed'),quote=F,col.names=F,sep="\t",row.names=F)


load('files/XOR.epiLCs.distinct.gene.interactions.Rdata')
genes.epiLCs<-unlist(lapply(as.character(XOR.epiLCs.distinct.gene.interactions[['gene']]),function(x) strsplit(x,",")))
genes.epiLCs<-gsub("\\..*","",genes.epiLCs)
epiLCs.distinct.coords<-unlist(lapply(as.character(XOR.epiLCs.distinct.gene.interactions[['XOR.epiLCs.coords']]),function(x) strsplit(x,",")))
epiLCs.distinct.distal.only<-unique(bedr.sort.region(unlist(lapply(epiLCs.distinct.coords,function(x) unlist(strsplit(x,"_"))[2]))))
epiLCs.distinct.df<-convert2bed(epiLCs.distinct.distal.only)
write.table(epiLCs.distinct.df,paste0('downstream/epiLCs.distinct.distal.bed'),quote=F,col.names=F,sep="\t",row.names=F)



