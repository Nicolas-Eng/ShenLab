bedtools intersect -a AD_full_list_srt.merged.bed -b EN.optimal.distal.res.atac.bed | sort -k1,1 -k2,2n - | bedtools merge -i - | wc -l

bedtools intersect -a AD_full_list_srt.merged.bed -b MG.optimal.distal.res.atac.bed | sort -k1,1 -k2,2n - | bedtools merge -i - | wc -l
