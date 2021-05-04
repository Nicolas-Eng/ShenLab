IDR.peaks<-read.table('ENCFF555XZZ_IDR_Peaks.bed',stringsAsFactors=FALSE,sep='\t')
IDR.peaks.thresh<-IDR.peaks[IDR.peaks$V7>2,][,1:3]

write.table(IDR.peaks.thresh,'ENCFF555XZZ_IDR_Peaks.0.01.bed',sep='\t',
	col.names = FALSE, row.names=FALSE,quote=FALSE)

library('rtracklayer')
library("GenomicRanges")

userRanges <- import.bed("ENCFF555XZZ_IDR_Peaks.0.01.bed")
resizeRanges <- resize(userRanges, width = 500)

export.bed(resizeRanges, "RightVentricle_Peaks_0.1.bed")

