mamba activate bioinfo

awk 'BEGIN{OFS="\t"}{print $0,NR}' mustache_Hansen_153k_mm39.bedpe > withID.bedpe

# anchor1: chr1, start1, end1, ID
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$7}' withID.bedpe > a1.bed
# anchor2: chr2, start2, end2, ID
awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$7}' withID.bedpe > a2.bed

liftOver -bedPlus=1 a1.bed /Volumes/UKJIN_SSD/Genomics_reference_common/Genome/liftOver_chainFiles/mm39ToMm10.over.chain.gz a1.mm10.bed a1.unmapped.bed
liftOver -bedPlus=1 a2.bed /Volumes/UKJIN_SSD/Genomics_reference_common/Genome/liftOver_chainFiles/mm39ToMm10.over.chain.gz a2.mm10.bed a2.unmapped.bed

# sort by ID (col4) so join works
sort -k4,4 a1.mm10.bed  > a1.sorted
sort -k4,4 a2.mm10.bed  > a2.sorted

# join on 4th column (ID), then spit out chr1,start1,end1,chr2,start2,end2
join -t$'\t' -1 4 -2 4 a1.sorted a2.sorted \
  | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5,$6,$7}' \
  > output.mm10.bedpe

sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n output.mm10.bedpe \
  > output.mm10.sorted.bedpe