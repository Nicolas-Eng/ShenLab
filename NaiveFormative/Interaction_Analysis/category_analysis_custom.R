library(bedr)
library(data.table)
library(plyr)

load('downstream_new_0.000001_test/trend.merged.df.Rdata')
load('files/XOR.interactions.gene.corrected.Rdata')


trend.type <- c('increase', 'decrease')

cell.types=c('mESCs','epiLCs')

trend.all.df <-list()
trend.XOR.interactions.all.df <-list()
genes.all.df <-list()

for (cell.type in cell.types) {
	cell.type.XOR<-XOR.interactions.gene.corrected[[cell.type]]
	for (trend2 in trend.type) {
		for(trend1 in trend.type) {
			trend.df<-trend.merged.df[trend.merged.df[['interaction.trend']]==trend1 & 
				trend.merged.df[['gene.trend']]==trend2 & 
				trend.merged.df[['DEseq']]=='yes',]
				genes<-trend.df[['gene']]
				print(head(trend.df,2))
				trend.XOR.interactions<-cell.type.XOR[which(cell.type.XOR[['gene']] %in% genes),]
				trend.distal.XOR.bed <- convert2bed(bedr.sort.region(unlist(trend.XOR.interactions[['distalcoord']])))
				oname<-paste0('downstream_new_0.000001_test/',cell.type,'.',trend1,'.',trend2,'.distal.bed')
#				write.table(trend.distal.XOR.bed,
#					file=oname,
#					quote=FALSE,row.names=FALSE,sep='\t',col.names=FALSE)
				print(oname)
				prefix<-paste0(cell.type,'.',trend1,'.',trend2)
				trend.all.df[[prefix]] = trend.df
				genes.all.df[[prefix]] = genes
				trend.XOR.interactions.all.df[[prefix]] = trend.XOR.interactions
				#save(trend.XOR.interactions,)
			}
		}
}



mESCs.di<-trend.all.df$mESCs.increase.increase
diff.MAPS<-as.numeric(mESCs.di[['epiLCs.fdr']])-as.numeric(mESCs.di[['mESCs.fdr']])
diff.RPKM<-as.numeric(mESCs.di[['epiLCs.RPKM']])-as.numeric(mESCs.di[['mESCs.RPKM']])

