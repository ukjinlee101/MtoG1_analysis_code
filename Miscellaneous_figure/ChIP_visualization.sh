# DRAWING METAPLOT & TORNADO PLOT ON UNION RAD21 PEAKS TO SHOW DEPLETION UPON dTAG TREATMENT
mamba activate macs2
chipDir="/Volumes/UKJIN_SSD/reference/ChIP"

# Merge peaks from ESC DMSO and ESC dTAG
cat ${chipDir}/RAD21_ESC-DMSO_unionPeaks.bed ${chipDir}/RAD21_ESC-dTAG_unionPeaks.bed \
| bedtools sort -i - \
| bedtools merge -d 200 -i - \
> "${chipDir}/RAD21_ESC_DMSOdTAG_unionPeaks.bed"


dataDir="/Volumes/UKJIN_SSD/data/deeptools"
figDir="/Volumes/UKJIN_SSD/figure/deeptools"
mkdir -p $dataDir
mkdir -p $figDir


name="RAD21_ESC-DMSOdTAG_unionPeaks"
computeMatrix reference-point \
	--referencePoint center \
	-R ${chipDir}/RAD21_ESC_DMSOdTAG_unionPeaks.bed \
	-S ${chipDir}/RAD21_ESC-DMSO_Rep1_spikeNorm_inputNorm.bw ${chipDir}/RAD21_ESC-DMSO_Rep2_spikeNorm_inputNorm.bw ${chipDir}/RAD21_ESC-dTAG_Rep1_spikeNorm_inputNorm.bw ${chipDir}/RAD21_ESC-dTAG_Rep2_spikeNorm_inputNorm.bw \
	--beforeRegionStartLength 2000 \
	--afterRegionStartLength 2000 \
	--binSize 25 \
	--missingDataAsZero \
	-p 8 \
	-o ${dataDir}/matrix_${name}.gz


plotHeatmap \
-m ${dataDir}/matrix_${name}.gz \
--heatmapHeight 0.5 \
--heatmapWidth 1.5 \
--colorMap Blues -max 1 \
-T $name --refPointLabel "" \
--samplesLabel "DMSO1" "DMSO2" "dTAG1" "dTAG2" \
--legendLocation none \
--yAxisLabel "spikeInputNorm" \
--xAxisLabel "peak center" \
--plotFileFormat svg \
-o ${figDir}/heatmap_${name}.svg


# Extracting RAD21 specific peaks
mamba activate bioinfo
chipDir="/Volumes/UKJIN_SSD/reference/ChIP"
BED1="${chipDir}/RAD21_ESC-DMSO_unionPeaks.bed"
BED2="${chipDir}/RAD21_ESC-dTAG_unionPeaks.bed"
BW="${chipDir}/RAD21_ESC-dTAG_Rep2_spikeNorm_inputNorm.bw"
TOPN="${4:-20}"
OUTDIR=$chipDir

BASENAME=$(basename "${BED2%.*}")

# 1) Peaks only in BED2, excluding any overlap with BED1
#    -v keeps A records that have no overlaps in B
bedtools intersect -v -a "$BED2" -b "$BED1" \
  > "$OUTDIR/${BASENAME}.only_in_dTAG.bed"

# 2) Convert bigWig to bedGraph
bigWigToBedGraph "$BW" "$OUTDIR/${BASENAME}.signal.bedgraph"

# 3) Map bigWig signal over peaks:
#    Append two columns to the BED: max and mean signal per peak
bedtools map \
  -a "$OUTDIR/${BASENAME}.only_in_dTAG.bed" \
  -b "$OUTDIR/${BASENAME}.signal.bedgraph" \
  -c 4 -o max,mean \
  > "$OUTDIR/${BASENAME}.only_in_dTAG.with_signal.bed"

# 4) Sort by max signal descending and take top N
#    The two appended columns are the last two fields; we sort on the penultimate field (max)
awk -v OFS='\t' '{print $0, (NF-1), NF}' "$OUTDIR/${BASENAME}.only_in_dTAG.with_signal.bed" | \
  sort -k$(awk 'END{print 0}' /dev/null) >/dev/null 2>&1 # placeholder to avoid shellcheck
# robust sort without awk gymnastics:
# extract max column index, then sort
MAXCOL=$(head -n1 "$OUTDIR/${BASENAME}.only_in_dTAG.with_signal.bed" | awk '{print NF-1}')
sort -k"${MAXCOL}","${MAXCOL}"gr \
  "$OUTDIR/${BASENAME}.only_in_dTAG.with_signal.bed" \
  | head -n "$TOPN" \
  > "$OUTDIR/${BASENAME}.top${TOPN}_by_max_signal.bed"

echo "Wrote:"
echo "  Peaks only in BED2:                $OUTDIR/${BASENAME}.only_in_bed2.bed"
echo "  Peaks with max,mean signal:        $OUTDIR/${BASENAME}.only_in_bed2.with_signal.bed"
echo "  Top ${TOPN} peaks by max signal:   $OUTDIR/${BASENAME}.top${TOPN}_by_max_signal.bed"