#cat WTC_*.narrowPeak > WTC_all_merged_q1e-4.shrt_peaks.bed 
#sort -k1,1 -k2,2n WTC_all_merged_q1e-4.shrt_peaks.bed | bedtools merge -i - | awk -vOFS='\t' '{$4="."; print}' - | awk -vOFS='\t' '$5=(FNR FS $5)' - | awk -vOFS='\t' '{$6="."; print}' - > WTC_all_merged_q1e-4.final.bed
#bgzip WTC_all_merged_q1e-4.final.bed
#tabix -p bed WTC_all_merged_q1e-4.final.bed.gz


sort -k1,1 -k2,2n WTC_0hr_merged.q1e-4.shrt_peaks.narrowPeak | bedtools merge -i - | awk -vOFS='\t' '{$4="."; print}' - | awk -vOFS='\t' '$5=(FNR FS $5)' - | awk -vOFS='\t' '{$6="."; print}' - > WTC_0hr_merged.q1e-4.final.bed
bgzip WTC_0hr_merged.q1e-4.final.bed
tabix -p bed WTC_0hr_merged.q1e-4.final.bed.gz

sort -k1,1 -k2,2n WTC_40hr_merged.q1e-4.shrt_peaks.narrowPeak | bedtools merge -i - | awk -vOFS='\t' '{$4="."; print}' - | awk -vOFS='\t' '$5=(FNR FS $5)' - | awk -vOFS='\t' '{$6="."; print}' - > WTC_40hr_merged.q1e-4.final.bed
bgzip WTC_40hr_merged.q1e-4.final.bed
tabix -p bed WTC_40hr_merged.q1e-4.final.bed.gz

sort -k1,1 -k2,2n WTC_75hr_merged.q1e-4.shrt_peaks.narrowPeak | bedtools merge -i - | awk -vOFS='\t' '{$4="."; print}' - | awk -vOFS='\t' '$5=(FNR FS $5)' - | awk -vOFS='\t' '{$6="."; print}' - > WTC_75hr_merged.q1e-4.final.bed
bgzip WTC_75hr_merged.q1e-4.final.bed
tabix -p bed WTC_75hr_merged.q1e-4.final.bed.gz