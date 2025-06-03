#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <sample name> <input.bw> <input.bed>"
  exit 1
fi

sample="$1"
bw="$2"
bed="$3"

bigWigToBedGraph ${bw} "${sample}.bg"
bedtools subtract -a "${sample}.bg" -b "${bed}" > "${sample}_masked.bg"
sort -k1,1V -k2,2n "${sample}_masked.bg" > "${sample}_masked_sorted.bg"
bedGraphToBigWig "${sample}_masked_sorted.bg" "/Volumes/UKJIN_SSD/MtoG1_analysis_code/reference/mm10/mm10.chrom.sizes" "${sample}_blacklistMasked.bw"

rm "${sample}.bg"
rm "${sample}_masked.bg"
rm "${sample}_masked_sorted.bg"