# From https://www.biostars.org/p/106249/
# Usage sh vcf_to_bed.sh intput.vcf > output.bed
bcftools norm -m- $1 | grep -v '^#' | awk -v OFS='\t' '{if(length($4) > length($5)) print $1,$2,$2+length($4)-1; else print $1,$2-1,$2+length($5)-1}' | bedtools sort -i - | bedtools merge -i -
