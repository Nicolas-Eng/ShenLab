cat EN.Interactions.x.ATAC.LHS.res.bed EN.Interactions.x.ATAC.RHS.res.bed | sort -k1,1 -k2,2n - | bedtools merge -i - > EN.Interactions.x.ATAC.all.res.srt.merged.bed

cat MG.Interactions.x.ATAC.LHS.res.bed MG.Interactions.x.ATAC.RHS.res.bed | sort -k1,1 -k2,2n - | bedtools merge -i - > MG.Interactions.x.ATAC.all.res.srt.merged.bed

bedtools intersect -a AD_full_list_srt.merged.bed -b EN.Interactions.x.ATAC.all.res.srt.merged.bed | sort -k1,1 -k2,2n - | bedtools merge -i - | wc -l

bedtools intersect -a AD_full_list_srt.merged.bed -b MG.Interactions.x.ATAC.all.res.srt.merged.bed | sort -k1,1 -k2,2n - | bedtools merge -i - | wc -l
