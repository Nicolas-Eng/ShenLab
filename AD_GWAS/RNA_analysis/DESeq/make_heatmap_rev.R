library(DESeq2)
library(tximport)
library(pheatmap)
library(RColorBrewer)
library(gdata)
library(biomaRt)
library(org.Hs.eg.db)
library(stringr)


MasterFile <- read.csv("DESeq_allReps_noEN.csv")
dir <- ("/shen/shenlabstore3/neng/20200309_MGEN/RNA_analysis/DESeq")
files <- file.path(dir, paste0(MasterFile$ID, ".merged.genes.results"))
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE,countsCol="expected_count")
txi.rsem$length[txi.rsem$length== 0] <- 1
ddsTxi <- DESeqDataSetFromTximport(txi.rsem, MasterFile, ~condition)


ddsTxi$condition<- as.factor(c("GFiMG Day2 rep1","GFiMG Day2 rep2",
	"GFiMG Day7 rep1","GFiMG Day7 rep2",
	"GFiMG Day14 rep1","GFiMG Day14 rep2",
	"GFiMG Day30 rep1","GFiMG Day30 rep2",
	"Ex vivo microglia (Gosselin et al. 2017) rep1","Ex vivo microglia (Gosselin et al. 2017) rep2","Ex vivo microglia (Gosselin et al. 2017) rep3",
	"In vitro Day7 microglia (Gosselin et al. 2017) rep1","In vitro Day7 microglia (Gosselin et al. 2017) rep2",
	"monocyte-derived macrophage (Heinz et al, 2018) rep1","monocyte-derived macrophage (Heinz et al, 2018) rep2",
	"monocyte-derived macrophage (Carlin et al, 2018) rep1","monocyte-derived macrophage (Carlin et al, 2018) rep2", 
	"WTC iPSC rep1", "WTC iPSC rep2"))

ddsTxi$ID<-as.factor(c("GFiMG Day2","GFiMG Day2",
	"GFiMG Day7","GFiMG Day7",
	"GFiMG Day14","GFiMG Day14",
	"GFiMG Day30","GFiMG Day30",
	"Ex vivo microglia (Gosselin et al. 2017)","Ex vivo microglia (Gosselin et al. 2017)","Ex vivo microglia (Gosselin et al. 2017)",
	"In vitro Day7 microglia (Gosselin et al. 2017)","In vitro Day7 microglia (Gosselin et al. 2017)",
	"monocyte-derived macrophage (Heinz et al, 2018)","monocyte-derived macrophage (Heinz et al, 2018)",
	"monocyte-derived macrophage (Carlin et al, 2018)","monocyte-derived macrophage (Carlin et al, 2018)", 
	"WTC iPSC", "WTC iPSC"))


mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
ens<-rownames(ddsTxi)
ensLookup <- gsub("\\.[0-9]*$", "", ens)
`%notin%` <- Negate(`%in%`)


MG_881genes<-read.table('genelist.txt')
ex<-as.character(MG_881genes[,1])
mart <- useEnsembl(biomart="ensembl", dataset = "hsapiens_gene_ensembl",GRCh=37)
gene2genomeEx <- getBM(values = ex,
  filters = "hgnc_symbol",
  mart = mart,
  attributes = c("hgnc_symbol", "ensembl_gene_id"))


genesFound<-gene2genomeEx[which(gene2genomeEx$ensembl_gene_id %in% ensLookup),]
genesNotFound<-ex[which(ex %notin% genesFound$hgnc_symbol)]

#Use Entrez to identify LOC IDs
genesLOC<-genesNotFound[str_detect(genesNotFound,'LOC')]

entrezLOC<-str_replace(genesLOC,"LOC","")
test <- getBM(values = entrezLOC,
  filters = "entrezgene_id",
  mart = mart,
  attributes = c("entrezgene_id","ensembl_gene_id"))
test$entrezgene_id<-paste('LOC',test$entrezgene_id,sep="")

genesFoundLOC<-test[which(test$ensembl_gene_id %in% ensLookup),]
LOCgene<-c('LOC102724323','LOC100505549')
Ensembleid<-c('ENSG00000231964','ENSG00000267040')

for (i in 1:length(LOCgene)) {
	genesFound<-rbind(genesFound,c(LOCgene[i],Ensembleid[i]))
}

colnames(genesFound)<-c("gene","ensembl_gene_id")
colnames(genesFoundLOC) <- c("gene","ensembl_gene_id")

genesFoundTot<-rbind(genesFound, genesFoundLOC)
###


####STEP 2#####
genesNotFoundPostLOC<-genesNotFound[!genesNotFound %in% genesLOC]

entrezmapping<-select(org.Hs.eg.db, 
   keys = genesNotFoundPostLOC,
   columns = c("ENTREZID", "SYMBOL"),
   keytype = "SYMBOL")

genesNotFoundPostLOCTest<-entrezmapping[which(!is.na(entrezmapping$ENTREZID)),]


noEntrezMappinghg38<-c('LINC01468','LINC01268','NCK1-AS1')
noEntrezMappinghg19<-c('RP11-346D6.6','RP1-249H1.4','RP11-85F14.5')
noEntrezEnsembl<-c('ENSG00000231131','ENSG00000227502','ENSG00000239213')

for (i in 1:length(noEntrezMappinghg19)) {
	genesFoundTot<-rbind(genesFoundTot,c(noEntrezMappinghg38[i],noEntrezEnsembl[i]))
}

test2 <- getBM(values = genesNotFoundPostLOCTest$ENTREZID,
  filters = "entrezgene_id",
  mart = mart,
  attributes = c("entrezgene_id","ensembl_gene_id"))

genesFound2<-test2[which(test2$ensembl_gene_id %in% ensLookup),]

for (i in 1:length(genesFound2$entrezgene_id)) {
	tmpID<-genesFound2$entrezgene_id[i]
	val<-which(genesNotFoundPostLOCTest$ENTREZID == tmpID)
	genesFound2$entrezgene_id[i]<-genesNotFoundPostLOCTest$SYMBOL[val]
}

colnames(genesFound2)<-c("gene","ensembl_gene_id")


missed_set<-genesNotFoundPostLOCTest[which(genesNotFoundPostLOCTest$SYMBOL %notin% genesFound2$gene),]
missedGenes2hg38<-c('FCGR1CP','LINC01410','LINC01907','LINC01909','LINC01480',
	'ADGRE2','SERPINB9P1','CATIP-AS1','LRP1-AS','STK26',
	'CEP295NL','ADGRG5','NRIR','TRPM2-AS',
	'CARD8-AS1','SUSD6','HACD4','RUBCNL',
	'UBR5-AS1','LINC00539','PYCARD-AS1','PDCD6IPP2','CXCL8',
	'P3H2','PLEKHM1P1','CFAP74','HLX-AS1','KLRA1P',
	'LMNTD2','MFSD4B','LINC01529','KANTR')
missedGenes2<-c('FCGR1C','RP11-262H14.1','AC098823.3','RP11-41O4.1','AC006129.4',
	'EMR2','RP11-420G6.4','AC021016.6','RP11-545N8.3','MST4',
	'DDC8','GPR114','AC017076.5','AP001065.2',
	'CTC-241F20.3','KIAA0247','PTPLAD2','KIAA0226L',
	'KB-431C1.4','LINC00422','C16orf98','RP11-578F21.12','IL8',
	'LEPREL1','PLEKHM1P','C1orf222','HLA-AS1','KLRA1P',
	'C11orf35','KIAA1919','AC002398.5','LINC01155')
ensembl_id2<-c('ENSG00000265531','ENSG00000238113','ENSG00000226125','ENSG00000266258','ENSG00000270164',
	'ENSG00000127507','ENSG00000230438','ENSG00000225062','ENSG00000259125','ENSG00000134602',
	'ENSG00000178404','ENSG00000159618','ENSG00000225964','ENSG00000230061',
	'ENSG00000268001', 'ENSG00000100647','ENSG00000188921','ENSG00000102445',
	'ENSG00000246263','ENSG00000224429','ENSG00000261359','ENSG00000261377','ENSG00000169429',
	'ENSG00000090530','ENSG00000214176','ENSG00000142609','ENSG00000257551','ENSG00000256667',
	'ENSG00000185522','ENSG00000173214','ENSG00000225872','ENSG00000232593')


for (i in 1:length(ensembl_id2)) {
	genesFound2<-rbind(genesFound2,c(missedGenes2hg38[i],ensembl_id2[i]))
}

genesFoundTot<-rbind(genesFoundTot,genesFound2)
print(genesFoundTot[which(genesFoundTot$ensembl_gene_id=='ENSG00000198019'),])

duplicates<-genesFoundTot[duplicated(genesFoundTot$ensembl_gene_id),]

number<-length(genesFoundTot[duplicated(genesFoundTot$ensembl_gene_id),][,1])


for (i in 1:number) {
	ensemble_check<-duplicates$ensembl_gene_id[i]
	print(genesFoundTot[which(genesFoundTot$ensembl_gene_id==ensemble_check),])
}



