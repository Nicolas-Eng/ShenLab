library(bedr)
library(data.table)
library(plyr)

args = commandArgs(trailingOnly=TRUE)
annotation.file<-args[1]
promoter.tmp.file<-args[2]
promoter.out.file<-args[3]


promoters<-read.table(gzfile(promoter.tmp.file,'rt'),header=F,stringsAsFactors=F)
gtf.annotation<-read.table(gzfile(annotation.file,'rt'),sep='\t',header=F,stringsAsFactors=F)

#Get Gene ID to Transcript ID 
#Get all rows that are transcripts
txids<-is.element(gtf.annotation$V3,"transcript")
#Get ninth column which has both gene_id and transcript_id
gtf.annotation.txids<-gtf.annotation[txids,]
#Split column
gtf.annotation.txids.split<-lapply(gtf.annotation.txids$V9,function(x) strsplit(x,'; '))
gene.ids<-unlist(lapply(gtf.annotation.txids.split,
	function(x) unlist(strsplit(unlist(x)[1]," "))[2]))
tx.ids<-unlist(lapply(gtf.annotation.txids.split,
	function(x) unlist(strsplit(unlist(x)[2]," "))[2]))

gene.tx.df<-data.frame(gene.ids,tx.ids)
colnames(promoters)<-c('chr','start','end','tx.ids','V5','strand')

promoters.gene.tx.ids.df<-merge(promoters,gene.tx.df,by="tx.ids")
reordered<-c('chr','start','end','gene.ids','tx.ids','V5','strand')
promoters.gene.tx.ids.df<-promoters.gene.tx.ids.df[reordered]

write.table(promoters.gene.tx.ids.df, promoter.out.file, sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)


