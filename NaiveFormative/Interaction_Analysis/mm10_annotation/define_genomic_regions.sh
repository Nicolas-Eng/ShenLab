gunzip -c gencode.vM25.annotation.gtf.gz | 
	awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' |
	sortBed |
	mergeBed -i - | gzip > gencode_vM25_exon_merged.bed.gz

gunzip -c gencode.vM25.annotation.gtf.gz |
	awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' |
	sortBed | subtractBed -a stdin -b gencode_vM25_exon_merged.bed.gz |
	gzip > gencode_vM25_intron.bed.gz

gunzip -c gencode.vM25.annotation.gtf.gz |
   awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' |
   sort -k1,1 -k2,2n - |
   complementBed -i stdin -g mm10.chrom.sizes |
   gzip > gencode_vM25_intergenic.bed.gz

check_utr.pl gencode.vM25.annotation.gtf.gz | gzip > transcript_utr_number.out.gz

perl promoter_mm10.pl gencode.vM25.annotation.gtf.gz 250 | gzip > promoter.tmp.bed.gz

Rscript promoter_addgenes.R gencode.vM25.annotation.gtf.gz promoter.tmp.bed.gz promoter.bed
sort -k1,1 -k2,2n promoter.bed | gzip - > promoter.bed.gz

rm promoter.tmp.bed.gz 