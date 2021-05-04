library(DiffBind)

analyze <- dba(sampleSheet = "/shen/shenlabstore3/neng/20200309_MGEN/DiffBind/MG_comparison_allReps.csv")
analyze <- dba.count(analyze, minOverlap=2,summits=100)
analyze <- dba.contrast(analyze, minMembers = 2)
analyze <- dba.analyze(analyze)

dba.report(analyze,file="MG_excneur_report_s100")

png('ATAC_heatmap.png', res=300, units="in", width=8, height=8)
dba.plotHeatmap(analyze)
dev.off()

png('ATAC_pca.png', res=300, units="in", width=8, height=8)
dba.plotPCA(analyze, label='ID')
dev.off()
