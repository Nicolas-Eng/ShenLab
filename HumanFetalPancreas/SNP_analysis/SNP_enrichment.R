library(bedr)
library(data.table)
library(plyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
library(qvalue)
library(annotables)
library(distributions3)
set.seed(1)
options(scipen=999)

cell.types<-c('RightVentricle','Alpha','Beta','Delta','Endocrine_Progenitor','Epsilon')
diseases<-c('T1D','T2D')

getBedrFormatted <- function(bed) {
	bedr <- paste0(bed[,1],':',bed[,2],'-',bed[,3])
	return(bedr)
}
 

#Format the RHS coordinates with bedr formatting chr:start-end
regroupAnnotationFile <- function(bed.inter) {
	tmpbed <- cbind(bed.inter[,2], bed.inter[,3], bed.inter[,4])
	#print(tmpbed)
	tmpbedr<-getBedrFormatted(tmpbed)
	bedr.inter <- cbind(bed.inter[[1]],tmpbedr)

	return(bedr.inter)
}

z.test = function(a, mu, var){
   zeta = (mean(a) - mu) / (sqrt(var / length(a)))
   return(zeta)
}

Z<-Normal(0,1)
SNP.enrichment.df<-data.frame(disease=character(),fg.val=numeric(),bg.stdev=numeric(),bg.mean=numeric(),
	z.score=numeric(),p.val=numeric(),significant=numeric(),cell.type=character(),stringsAsFactors=FALSE)


disease.SNPs.list<-list()

for (disease in diseases) {
	disease.SNPs<-read.table(paste0('/shen/shenlabstore3/neng/20200509_scATAC_sneddon/SNP_analysis/features/',disease,'.complete.snps.features.bed'),header=F, stringsAsFactors=F)
	disease.SNPs$V1<-paste0('chr',disease.SNPs$V1)
	disease.SNPs.list[[disease]]<-disease.SNPs
}

save(disease.SNPs.list,file='saved/disease.SNPs.list.Rdata')

bedr.SNPs.list<-list()
bedr.Peaks.list<-list()
for (disease in diseases) {
	disease.SNPs<-read.table(paste0('/shen/shenlabstore3/neng/20200509_scATAC_sneddon/SNP_analysis/features/',disease,'.complete.snps.features.bed'),header=F, stringsAsFactors=F)
	disease.SNPs$V1<-paste0('chr',disease.SNPs$V1)

	for (cell.type in cell.types) {
		cell.type.peaks<-read.table(paste0('/shen/shenlabstore3/neng/20200509_scATAC_sneddon/SNP_analysis/12wpc_endocrine_peaks/',cell.type,'/',cell.type,'_Peaks_0.1.bed'), header=F, stringsAsFactors=F)
		bedr.peaks.SNPs.intersect<-bedr(input=list(a=bedr.sort.region(getBedrFormatted(disease.SNPs)),
			b=bedr.sort.region(getBedrFormatted(cell.type.peaks))),method="intersect",params="-u")
		bedr.peaks.SNPs.regions<-bedr(input=list(b=bedr.sort.region(getBedrFormatted(disease.SNPs)),
			a=bedr.sort.region(getBedrFormatted(cell.type.peaks))),method="intersect",params="-u")


		fg.val<-length(bedr.peaks.SNPs.intersect)
		bedr.SNPs.list[[disease]][[cell.type]]<-bedr.peaks.SNPs.intersect
		bedr.Peaks.list[[disease]][[cell.type]]<-bedr.peaks.SNPs.regions

		print(fg.val)
		bg.list<-c()
		for (i in 2:101) {
			ctrl.type.peaks<-read.table(paste0('/shen/shenlabstore3/neng/20200509_scATAC_sneddon/SNP_analysis/12wpc_endocrine_control/',cell.type,'_Peaks_ctrl_s',i,'.bed'),header=F,stringsAsFactors=F)
			bedr.ctrl.peaks.SNPs.intersect<-bedr(input=list(a=bedr.sort.region(getBedrFormatted(disease.SNPs)),
				b=bedr.sort.region(getBedrFormatted(ctrl.type.peaks))),method="intersect",params="-u")
			print(length(bedr.ctrl.peaks.SNPs.intersect))
			bg.list<-c(bg.list,length(bedr.ctrl.peaks.SNPs.intersect))
			print(paste(i,'ctrl number'))
		}
		bg.stdev<-sd(bg.list)/ sqrt(length(bg.list))
		bg.mean<-mean(bg.list)
		z.score = (fg.val - bg.mean) / (bg.stdev)
		p.val= pnorm(z.score, lower.tail=FALSE)
		print(p.val)



		significant=''
		if (p.val <= 0.05) {
			significant='yes'
		} else {
			significant='no'
		}

		new.row<-data.frame(disease=as.character(disease), fg.val=fg.val, bg.stdev=bg.stdev, 
			bg.mean=bg.mean,z.score=z.score,p.val=p.val, 
			significant=as.character(significant),cell.type=as.character(cell.type))
		print(new.row)
		SNP.enrichment.df<-rbind(SNP.enrichment.df,new.row)
		print(SNP.enrichment.df)
	}

}

save(bedr.SNPs.list,file='saved/bedr.SNPs.list.Rdata')
save(bedr.Peaks.list,file='saved/bedr.Peaks.list.Rdata')
save(SNP.enrichment.df,file='saved/SNP.enrichment.df.Rdata')
load('saved/bedr.Peaks.list.Rdata')
load('saved/SNP.enrichment.df.Rdata')
load('saved/bedr.SNPs.list.Rdata')
library(ggplot2)
library(ggsignif)
library(VennDiagram)

SNP.enrichment.ggplot2<-data.frame(disease=rep(SNP.enrichment.df$disease,2),
	celltype=rep(SNP.enrichment.df$cell.type,2),
	value=as.numeric(c(SNP.enrichment.df$fg.val,SNP.enrichment.df$bg.mean)),
	type=c(rep('ATAC',length(SNP.enrichment.df$disease)),rep('controls',length(SNP.enrichment.df$disease))),
	standarderror=as.numeric(c(rep(0,length(SNP.enrichment.df$disease)),SNP.enrichment.df$bg.stdev)))


for (disease in diseases) {
	colors=''
	if (disease=='T1D') {
		colors=c('#5DADE2','#D6EAF8')
		y.pos=c(93, 119, 140, 131, 106, 82)
		p.val=c('NS',rep('<2.2E-16',5))
		ylim=c(0,150)
	} else {
		colors=c('#239B56','#D5F5E3')
		y.pos=c(178, 305, 334, 309, 290, 229)
		p.val=c('NS',rep('<2.2E-16',5))
		ylim=c(0,350)
	}
	SNP.enrichment.disease<-SNP.enrichment.ggplot2[SNP.enrichment.ggplot2$disease==disease,]
	f <- ggplot(SNP.enrichment.disease, aes(x=celltype, y=value, fill=type))
	f <- f + geom_bar(position="dodge", size=1,stat="identity", colour='black')
	f <- f + geom_errorbar(aes(ymin=value-standarderror, ymax=value+standarderror), width=.2,position=position_dodge(.9))
	f <- f + theme_classic()
	f <- f + scale_fill_manual(values=colors)
	f <- f + ylab(paste0('Number of ',disease,' SNPs')) + xlab('')
	f <- f + scale_y_continuous(expand = c(0, 0), limits = ylim)  + scale_x_discrete(labels= c('Right Ventricle', 'Alpha', 'Beta', 'Delta','Endocrine Progenitor', 'Epsilon'))
	f <- f + geom_signif(y_position=y.pos, xmin=c(0.7, 1.7, 2.7, 3.7, 4.7, 5.7), xmax=c(1.3, 2.3, 3.3, 4.3, 5.3, 6.3),
              annotation=p.val) 




	venn.file.name<-''

	if (disease=='T1D') {
		ggsave('T1D.SNP.Enrichment.pdf',width=9,height=9)
		venn.file.name<-'T1D_SNP_Venn.png'

	} else {
		ggsave('T2D.SNP.Enrichment.pdf',width=9,height=9)
		venn.file.name<-'T2D_SNP_Venn.png'
	}


	venn.diagram(x = list(Alpha=bedr.SNPs.list[[disease]]$Alpha,
		Beta=bedr.SNPs.list[[disease]]$Beta,
		Delta=bedr.SNPs.list[[disease]]$Delta,
		EP=bedr.SNPs.list[[disease]]$Endocrine_Progenitor,
		Epsilon=bedr.SNPs.list[[disease]]$Epsilon),
		filename = venn.file.name,
		col = "black",
		fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
		alpha = 0.50,
		cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
		 1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
		cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
		cat.cex = 1.5,
		cat.fontface = "bold",
		margin = 0.05,
		height=5000,
		width=5000
		);


}


### check if bedr.intersected peaks are in the original bed file ###
for(disease in diseases) {
	for (cell.type in cell.types) {
		cell.type.peaks<-getBedrFormatted(read.table(paste0('/shen/shenlabstore3/neng/20200509_scATAC_sneddon/SNP_analysis/12wpc_endocrine_peaks/',cell.type,'/',cell.type,'_Peaks_0.1.bed'), header=F, stringsAsFactors=F))
		bedr.peaks<-bedr.Peaks.list[[disease]][[cell.type]]
		if(bedr.peaks %in% cell.type.peaks) {
			print(bedr.peaks %in% cell.type.peaks)
		}
	}
}





