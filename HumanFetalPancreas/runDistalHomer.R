library(rtracklayer)
FEVhigh<-readRDS('Save-ProjEndo2/PeakCalls/FEV.High-reproduciblePeaks.gr.rds')
FEVhighNC<-FEVhigh[!FEVhigh$peakType=="Promoter" & !FEVhigh$peakType=="Exonic"]
export(FEVhighNC, format='bed',"HOMER_analysis/FEV.noncoding.bed")


Beta<-readRDS('Save-ProjEndo2/PeakCalls/Beta-reproduciblePeaks.gr.rds')
BetaNC<-Beta[!Beta$peakType=="Promoter" & !Beta$peakType=="Exonic"]
export(BetaNC, format='bed',"HOMER_analysis/Beta.noncoding.bed")

Common_Progenitor<-readRDS('Save-ProjEndo2/PeakCalls/Common_Progenitor-reproduciblePeaks.gr.rds')
Common_ProgenitorNC<-Common_Progenitor[!Common_Progenitor$peakType=="Promoter" & !Common_Progenitor$peakType=="Exonic"]
export(Common_ProgenitorNC, format='bed',"HOMER_analysis/Common_Progenitor.noncoding.bed")