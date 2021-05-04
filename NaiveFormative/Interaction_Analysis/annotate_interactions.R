library(bedr)
library(data.table)
library(plyr)


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


annotateXORinteractions <- function(XOR.interactions, anchors) {
	XOR.LHS.interactions<-getBedrFormatted(XOR.interactions[1:3])
	XOR.RHS.interactions<-getBedrFormatted(XOR.interactions[4:6])

	XOR.reference.annotated<-cbind(XOR.LHS.interactions,XOR.RHS.interactions,XOR.interactions[9])

	#Sort Files
	XOR.LHS.interactions.sort <- bedr.sort.region(XOR.LHS.interactions)
	XOR.RHS.interactions.sort <- bedr.sort.region(XOR.RHS.interactions)
	anchors.sort<-bedr.sort.region(anchors.annotated.df[[1]])

	bed.XOR.LHS.intersected <- bedr( input = list(a=XOR.LHS.interactions.sort, b=anchors.sort), method = "intersect", params = "-wo")
	bedr.XOR.LHS.intersected <- regroupAnnotationFile(bed.XOR.LHS.intersected)

	bed.XOR.RHS.intersected <- bedr( input = list(a=XOR.RHS.interactions.sort, b=anchors.sort), method = "intersect", params = "-wo")
	bedr.XOR.RHS.intersected <- regroupAnnotationFile(bed.XOR.RHS.intersected)

	colnames(bedr.XOR.LHS.intersected)<-c('XOR.LHS.interactions','anchors')
	colnames(bedr.XOR.RHS.intersected)<-c('XOR.RHS.interactions','anchors')


	bedr.XOR.RHS.intersected.annotated<-join(data.frame(bedr.XOR.RHS.intersected),data.frame(anchors.annotated.df),"anchors")
	bedr.XOR.LHS.intersected.annotated<-join(data.frame(bedr.XOR.LHS.intersected),data.frame(anchors.annotated.df),"anchors")


	bedr.XOR.RHS.collapsed<-collapseAnnotations(bedr.XOR.RHS.intersected.annotated,"XOR.RHS.interactions",c("gene"))
	bedr.XOR.LHS.collapsed<-collapseAnnotations(bedr.XOR.LHS.intersected.annotated,"XOR.LHS.interactions",c("gene"))
	colnames(bedr.XOR.LHS.collapsed)<-c('XOR.LHS.interactions','LHS.anchor.genes')
	colnames(bedr.XOR.RHS.collapsed)<-c('XOR.RHS.interactions','RHS.anchor.genes')

	tmpRHS<-na.omit(join(XOR.reference.annotated,bedr.XOR.RHS.collapsed,'XOR.RHS.interactions'))
	tmpLHS<-na.omit(join(XOR.reference.annotated,bedr.XOR.LHS.collapsed,'XOR.LHS.interactions'))

	#Write format so that it is: anchor, distal, fdr
	colnames(tmpRHS)<-c("distalcoords","anchorcoords","fdr","gene")
	revtmpRHS<-tmpRHS[c("anchorcoords","distalcoords","fdr","gene")]
	colnames(tmpLHS)<-c("anchorcoords","distalcoords","fdr","gene")

	bedr.XOR.annotated.ordered<-rbind(revtmpRHS,tmpLHS)
	bedr.XOR.annotated.ordered[[1]]<-as.character(bedr.XOR.annotated.ordered[[1]])
	bedr.XOR.annotated.ordered[[2]]<-as.character(bedr.XOR.annotated.ordered[[2]])
	return(bedr.XOR.annotated.ordered)
}

regroupInteractionAnnotation<-function(interaction.file,cell.type) {
	XOR.interactions<-interaction.file[[cell.type]]
	XOR.coords<-paste0(XOR.interactions[["anchorcoords"]],'_',XOR.interactions[["distalcoords"]])
	XOR.combined<-cbind(XOR.coords,XOR.interactions[3:4],stringsAsFactors=FALSE)
	colnames(XOR.combined)<-c(paste0('XOR.',cell.type,'.coords'),paste0(cell.type,'.fdr'),'gene')
	return(XOR.combined)
}

load('files/anchors.annotated.Rdata')
cell.types<-c('mESCs','epiLCs')

XOR.interactions.annotated<-list()
for (cell.type in cell.types) {
	XOR.df<-read.table(paste0(cell.type,'_H3K4me3_trim_final_merged.5k.2.peaks.xor'),sep="\t", header=T, stringsAsFactors=F)
	XOR.cell.type.annotated<-annotateXORinteractions(XOR.df,anchors.annotated.df)
	XOR.interactions.annotated[[cell.type]]<-XOR.cell.type.annotated
}

save(XOR.interactions.annotated,file='files/XOR.interactions.annotated.Rdata')
load('files/XOR.interactions.annotated.Rdata')

XOR.interactions.annotated.fixed<-list()
for (cell.type in cell.types) {
	XOR.cell.type.combined<-regroupInteractionAnnotation(XOR.interactions.annotated,cell.type)
	XOR.interactions.annotated.fixed[[cell.type]]<-XOR.cell.type.combined
}


XOR.combined<-dplyr::bind_rows(XOR.interactions.annotated.fixed[['mESCs']],XOR.interactions.annotated.fixed[['epiLCs']])
XOR.combined.gene.collapse<-collapseAnnotations(XOR.combined,"gene",c("XOR.epiLCs.coords","XOR.mESCs.coords","mESCs.fdr","epiLCs.fdr"))

#Get genes interacted with only mESCs or epiLCs
XOR.epiLCs.distinct.gene.interactions<-XOR.combined.gene.collapse[XOR.combined.gene.collapse[["XOR.mESCs.coords"]]=="NA",]
XOR.mESCs.distinct.gene.interactions<-XOR.combined.gene.collapse[XOR.combined.gene.collapse[["XOR.epiLCs.coords"]]=="NA",]

save(XOR.epiLCs.distinct.gene.interactions,file='files/XOR.epiLCs.distinct.gene.interactions.Rdata')
save(XOR.mESCs.distinct.gene.interactions,file='files/XOR.mESCs.distinct.gene.interactions.Rdata')

#Get genes interacted with both cell-types
XOR.shared.gene.interactions<-XOR.combined.gene.collapse[XOR.combined.gene.collapse[["XOR.mESCs.coords"]]!="NA" 
	& XOR.combined.gene.collapse[["XOR.epiLCs.coords"]]!="NA",]

#Turn each coordinate pair into a item of a list for comparison
XOR.shared.gene.interactions[["XOR.mESCs.coords"]]<-lapply(XOR.shared.gene.interactions[["XOR.mESCs.coords"]],
	function(x) unlist(strsplit(x,","))[unlist(strsplit(x,",")) != "NA"])
XOR.shared.gene.interactions[["XOR.epiLCs.coords"]]<-lapply(XOR.shared.gene.interactions[["XOR.epiLCs.coords"]],
	function(x) unlist(strsplit(x,","))[unlist(strsplit(x,",")) != "NA"])
XOR.shared.gene.interactions[["mESCs.fdr"]]<-lapply(XOR.shared.gene.interactions[["mESCs.fdr"]],
	function(x) unlist(strsplit(x,","))[unlist(strsplit(x,",")) != "NA"])
XOR.shared.gene.interactions[["epiLCs.fdr"]]<-lapply(XOR.shared.gene.interactions[["epiLCs.fdr"]],
	function(x) unlist(strsplit(x,","))[unlist(strsplit(x,",")) != "NA"])

#Classify interactions for each gene as shared, mESCs, or epiLCs only
XOR.classify.shared.gene.interactions<-data.frame()

for (i in 1:length(XOR.shared.gene.interactions[,1])){
	XOR.gene<-unlist(XOR.shared.gene.interactions[i,"gene"])
	XOR.mESCs.interaction<-unlist(XOR.shared.gene.interactions[i,"XOR.mESCs.coords"])
	XOR.epiLCs.interaction<-unlist(XOR.shared.gene.interactions[i,"XOR.epiLCs.coords"])
	XOR.shared.interaction<-intersect(XOR.mESCs.interaction,XOR.epiLCs.interaction)
	XOR.mESCs.only.interaction<-setdiff(XOR.mESCs.interaction,XOR.epiLCs.interaction)
	XOR.epiLCs.only.interaction<-setdiff(XOR.epiLCs.interaction,XOR.mESCs.interaction)
	XOR.shared.col<-gsub(" ","",toString(XOR.shared.interaction))
	XOR.mESCs.only.col<-gsub(" ","",toString(XOR.mESCs.only.interaction))
	XOR.epiLCs.only.col<-gsub(" ","",toString(XOR.epiLCs.only.interaction))
	XOR.separated<-cbind(XOR.gene,XOR.shared.col,XOR.mESCs.only.col,XOR.epiLCs.only.col)
	if (empty(XOR.classify.shared.gene.interactions)){
		XOR.classify.shared.gene.interactions<-XOR.separated
	} else {
		XOR.classify.shared.gene.interactions<-rbind(XOR.classify.shared.gene.interactions,XOR.separated)
	}
}


XOR.csgi.df<-data.frame(XOR.classify.shared.gene.interactions)
colnames(XOR.csgi.df)<-c('gene','shared','mESCs.only','epiLCs.only')

XOR.shared.no.switch.gene.interactions<-XOR.csgi.df[(XOR.csgi.df['mESCs.only']=="" & XOR.csgi.df['epiLCs.only']==""),]
XOR.shared.switch.gene.interactions<-XOR.csgi.df[!(XOR.csgi.df['mESCs.only']=="" & XOR.csgi.df['epiLCs.only']==""),]

save(XOR.shared.no.switch.gene.interactions,file='files/XOR.shared.no.switch.gene.interactions.Rdata')
save(XOR.shared.switch.gene.interactions,file='files/XOR.shared.switch.gene.interactions.Rdata')



