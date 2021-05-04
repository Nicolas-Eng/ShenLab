cat ADSNP_imputeTOPMADexpand_r2_0.8_hg19.bed hg19_17250AD.bed > AD_full_list.bed
sort -k1,1 -k2,2n AD_full_list.bed > AD_full_list_srt.bed
#mergeBed -i AD_full_list_srt.bed -c 4,4 -o collapse, count > AD_full_list_srt.merged.bed
mergeBed -i AD_full_list_srt.bed -c 4 -o collapse > AD_full_list_srt.merged.bed
bedtools intersect -a AD_full_list_srt.merged.bed -b hg19.baitmap -wo > AD_full_SNPs_x_baitmap.bed
cut -f5-7 AD_full_SNPs_x_baitmap.bed | bedtools merge -i - -c 1 -o count > AD_full_SNPs_x_baitmap_coords_count.bed

bedtools intersect -a AD_full_SNPs_x_baitmap_coords_count.bed -b ../distalATAC/MG.optimal.distal.atac.bed > AD_SNPs_baitmap_coords_count_x_MG_distal.atac.bed
bedtools intersect -a AD_full_SNPs_x_baitmap_coords_count.bed -b ../distalATAC/EN.optimal.distal.atac.bed > AD_SNPs_baitmap_coords_count_x_EN_distal.atac.bed