library(bedr)
library(data.table)
library(plyr)
anchors<-read.table('anchors.bed',sep="\t", header=F, stringsAsFactors=F)
promoters<-read.table(gzfile('../mm10_annotation/promoter.bed.gz','rt'),header=F,stringsAsFactors=F)


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


#Format for bedr Manipulation
anchors_bedr <- getBedrFormatted(anchors)
promoters_bedr <- getBedrFormatted(promoters)

#Add Back Transcript Annotation
promoters_annotated<-cbind(promoters_bedr,promoters[,4],promoters[,5])
colnames(promoters_annotated)<-c('coords','gene','transcript')

#Sort Files
anchors.sort <- bedr.sort.region(anchors_bedr)
promoters.sort <- bedr.sort.region(promoters_bedr)

#Intersect anchors with promoter bed file
abed.intersected <- bedr( input = list(a=anchors.sort, b=promoters.sort), method = "intersect", params = "-wo")
abedr.intersected <- regroupAnnotationFile(abed.intersected)
colnames(abedr.intersected) <- c('anchors','coords')

#Merge Data.Frames by Coordinates
abedr.intersected.annotated<-join(data.frame(abedr.intersected),data.frame(promoters_annotated),'coords')

anchors.annotated.df<-collapseAnnotations(abedr.intersected.annotated,"anchors",c("gene","transcript"))
anchors.gene.annotated.df<-collapseAnnotations(abedr.intersected.annotated,"gene",c("anchors","transcript"))

save(anchors.gene.annotated.df,file='files/anchors.gene.annotated.Rdata')
save(anchors.annotated.df,file='files/anchors.annotated.Rdata')
#load('files/anchors.annotated.Rdata')
