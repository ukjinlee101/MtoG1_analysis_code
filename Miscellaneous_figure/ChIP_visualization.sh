############################################################################################
# REVIEWED FOR PUBLICATION START
# DRAWING METAPLOT & TORNADO PLOT ON UNION RAD21 PEAKS TO SHOW DEPLETION UPON dTAG TREATMENT
mamba activate macs2
chipDir="/Volumes/UKJIN_SSD/reference/ChIP/RAD21-HA"

# Merge peaks from ESC DMSO and ESC dTAG
peak1="${chipDir}/RAD21_ESC-DMSO_unionPeaks.bed"
peak2="${chipDir}/RAD21_ESC-dTAG_unionPeaks.bed"
out="${chipDir}/UnionPeaks_RAD21_ESC-DMSO_and_dTAG.bed"
if [[ -f ${peak1} && -f ${peak2} ]]; then
    cat ${peak1} ${peak2} \
    | bedtools sort -i - \
    | bedtools merge -d 200 -i - \
    > "$out"
else
    echo "Missing input files"
fi

# Filtering union peaks that are present in each condition
out1="${chipDir}/UnionPeaks_RAD21_ESC-DMSO_and_dTAG__ESC-DMSO.bed"
bedtools intersect -wa -a ${out} -b ${peak1} > ${out1}
out2="${chipDir}/UnionPeaks_RAD21_ESC-DMSO_and_dTAG__ESC-dTAG.bed"
bedtools intersect -wa -a ${out} -b ${peak2} > ${out2}
out3="${chipDir}/UnionPeaks_RAD21_ESC-DMSO_and_dTAG__common.bed"
bedtools intersect -wa -a ${out1} -b ${out2} > ${out3}

############################################################################################
# Drawing tornado plot
dataDir="/Volumes/UKJIN_SSD/data/deeptools"
figDir="/Volumes/UKJIN_SSD/figure/deeptools"
mkdir -p $dataDir
mkdir -p $figDir


name="UnionPeaks_RAD21_ESC-DMSO_and_dTAG"
computeMatrix reference-point \
	--referencePoint center \
	-R ${chipDir}/UnionPeaks_RAD21_ESC-DMSO_and_dTAG.bed \
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



name="UnionPeaks_RAD21_ESC-DMSO_and_dTAG__common"
computeMatrix reference-point \
    --referencePoint center \
    -R ${chipDir}/UnionPeaks_RAD21_ESC-DMSO_and_dTAG__common.bed \
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


############################################################################################
# Merge peaks from HA and RAD21 for ESC and EpiLC
peak1="${chipDir}/RAD21_ESC-DMSO_unionPeaks.bed"
peak2="${chipDir}/HA_ESC-DMSO_unionPeaks.bed"
out="${chipDir}/UnionPeaks_RAD21-HA_ESC-DMSO.bed"
if [[ -f ${peak1} && -f ${peak2} ]]; then
    cat ${peak1} ${peak2} \
    | bedtools sort -i - \
    | bedtools merge -d 200 -i - \
    > "$out"
else
    echo "Missing input files"
fi

peak1="${chipDir}/RAD21_EpiLC_unionPeaks.bed"
peak2="${chipDir}/HA_EpiLC_unionPeaks.bed"
out="${chipDir}/UnionPeaks_RAD21-HA_EpiLC.bed"
if [[ -f ${peak1} && -f ${peak2} ]]; then
    cat ${peak1} ${peak2} \
    | bedtools sort -i - \
    | bedtools merge -d 200 -i - \
    > "$out"
else
    echo "Missing input files"
fi

peak1="${chipDir}/UnionPeaks_RAD21-HA_ESC-DMSO.bed"
peak2="${chipDir}/UnionPeaks_RAD21-HA_EpiLC.bed"
out="${chipDir}/UnionPeaks_RAD21-HA_ESC-DMSO_and_EpiLC.bed"
if [[ -f ${peak1} && -f ${peak2} ]]; then
    cat ${peak1} ${peak2} \
    | bedtools sort -i - \
    | bedtools merge -d 200 -i - \
    > "$out"
else
    echo "Missing input files"
fi

# Filtering union peaks that are present in each condition
out="${chipDir}/UnionPeaks_RAD21-HA_ESC-DMSO_and_EpiLC.bed"
peak1="${chipDir}/UnionPeaks_RAD21-HA_ESC-DMSO.bed"
peak2="${chipDir}/UnionPeaks_RAD21-HA_EpiLC.bed"
out1="${chipDir}/UnionPeaks_RAD21-HA_ESC-DMSO_and_EpiLC__ESC-DMSO.bed"
bedtools intersect -wa -a ${out} -b ${peak1} > ${out1}
out2="${chipDir}/UnionPeaks_RAD21-HA_ESC-DMSO_and_EpiLC__EpiLC.bed"
bedtools intersect -wa -a ${out} -b ${peak2} > ${out2}
out3="${chipDir}/UnionPeaks_RAD21-HA_ESC-DMSO_and_EpiLC__common.bed"
bedtools intersect -wa -a ${out1} -b ${out2} > ${out3}


############################################################################################

dataDir="/Volumes/UKJIN_SSD/data/deeptools"
figDir="/Volumes/UKJIN_SSD/figure/deeptools"
mkdir -p $dataDir
mkdir -p $figDir


name="UnionPeaks_RAD21-HA_ESC-DMSO_and_EpiLC"
computeMatrix reference-point \
  --referencePoint center \
  -R ${chipDir}/UnionPeaks_RAD21-HA_ESC-DMSO_and_EpiLC.bed \
  -S ${chipDir}/RAD21_ESC-DMSO_Rep1_spikeNorm_inputNorm.bw ${chipDir}/RAD21_ESC-DMSO_Rep2_spikeNorm_inputNorm.bw ${chipDir}/HA_ESC-DMSO_Rep1_spikeNorm_inputNorm.bw ${chipDir}/HA_ESC-DMSO_Rep2_spikeNorm_inputNorm.bw ${chipDir}/RAD21_EpiLC_Rep1_spikeNorm_inputNorm.bw ${chipDir}/RAD21_EpiLC_Rep2_spikeNorm_inputNorm.bw ${chipDir}/HA_EpiLC_Rep1_spikeNorm_inputNorm.bw ${chipDir}/HA_EpiLC_Rep2_spikeNorm_inputNorm.bw \
  --beforeRegionStartLength 2000 \
  --afterRegionStartLength 2000 \
  --binSize 25 \
  --missingDataAsZero \
  -p 8 \
  -o ${dataDir}/matrix_${name}.gz


plotHeatmap \
-m ${dataDir}/matrix_${name}.gz \
--heatmapHeight 0.5 \
--heatmapWidth 3 \
--colorMap Blues Blues Greens Greens Blues Blues Greens Greens -max 1 \
-T $name --refPointLabel "" \
--samplesLabel "ESC1" "ESC2" "ESC3" "ESC4" "EpiLC1" "EpiLC2" "EpiLC3" "EpiLC4" \
--legendLocation none \
--yAxisLabel "spikeInputNorm" \
--xAxisLabel "peak center" \
--plotFileFormat svg \
-o ${figDir}/heatmap_${name}.svg


############################################################################################

dataDir="/Volumes/UKJIN_SSD/data/deeptools"
figDir="/Volumes/UKJIN_SSD/figure/deeptools"
mkdir -p $dataDir
mkdir -p $figDir


name="UnionPeaks_RAD21-HA_ESC-DMSO_and_EpiLC"
computeMatrix reference-point \
  --referencePoint center \
  -R ${chipDir}/UnionPeaks_RAD21-HA_ESC-DMSO_and_EpiLC.bed \
  -S ${chipDir}/RAD21_ESC-DMSO_Rep1_spikeNorm_inputNorm.bw ${chipDir}/RAD21_ESC-DMSO_Rep2_spikeNorm_inputNorm.bw ${chipDir}/HA_ESC-DMSO_Rep1_spikeNorm_inputNorm.bw ${chipDir}/HA_ESC-DMSO_Rep2_spikeNorm_inputNorm.bw ${chipDir}/RAD21_EpiLC_Rep1_spikeNorm_inputNorm.bw ${chipDir}/RAD21_EpiLC_Rep2_spikeNorm_inputNorm.bw ${chipDir}/HA_EpiLC_Rep1_spikeNorm_inputNorm.bw ${chipDir}/HA_EpiLC_Rep2_spikeNorm_inputNorm.bw \
  --beforeRegionStartLength 2000 \
  --afterRegionStartLength 2000 \
  --binSize 25 \
  --missingDataAsZero \
  -p 8 \
  -o ${dataDir}/matrix_${name}.gz


plotHeatmap \
-m ${dataDir}/matrix_${name}.gz \
--heatmapHeight 0.5 \
--heatmapWidth 3 \
--colorMap Blues Blues Greens Greens Blues Blues Greens Greens -max 1 \
-T $name --refPointLabel "" \
--samplesLabel "ESC1" "ESC2" "ESC3" "ESC4" "EpiLC1" "EpiLC2" "EpiLC3" "EpiLC4" \
--legendLocation none \
--yAxisLabel "spikeInputNorm" \
--xAxisLabel "peak center" \
--plotFileFormat svg \
-o ${figDir}/heatmap_${name}.svg



# DRAWING METAPLOT TO SHOW THE REPRODUCIBILITY ACROSS REPLICATES
mamba activate deeptools

workDir="/Volumes/UKJIN_SSD/ChIP_250609_18886_H3K4me3K27ac/pipeline/result"
outDir="${workDir}/../metaplot"
mkdir -p $outDir

# Diff plot
computeMatrix reference-point \
--referencePoint center \
-R H3K4me3_diff_specific_ESC.bed H3K4me3_diff_specific_EpiLC.bed \
-S ${workDir}/bigwig/H3K4me3_ESC_Rep1_depthNorm.bw ${workDir}/bigwig/H3K4me3_ESC_Rep2_depthNorm.bw ${workDir}/bigwig/H3K4me3_EpiLC_Rep1_depthNorm.bw ${workDir}/bigwig/H3K4me3_EpiLC_Rep2_depthNorm.bw \
--beforeRegionStartLength 2000 \
--afterRegionStartLength 2000 \
--binSize 25 \
--missingDataAsZero \
-p 8 \
-o ${outDir}/matrix_diff_H3K4me3_onlyChange.gz

plotHeatmap \
-m ${outDir}/matrix_diff_H3K4me3_onlyChange.gz \
--heatmapHeight 0.5 \
--heatmapWidth 3 \
--colorMap Blues Blues Greens Greens \
-T "" --refPointLabel "" \
--samplesLabel "ESC1" "ESC2" "EpiLC1" "EpiLC2" \
--legendLocation none \
--yAxisLabel "RPM" \
--xAxisLabel "peak center" \
--regionsLabel "ESC" "EpiLC" \
-o ${outDir}/heatmap_diff_H3K4me3_onlyChange.svg


computeMatrix reference-point \
--referencePoint center \
-R H3K27ac_diff_specific_ESC.bed H3K27ac_diff_specific_EpiLC.bed \
-S ${workDir}/bigwig/H3K27ac_ESC_Rep1_depthNorm.bw ${workDir}/bigwig/H3K27ac_ESC_Rep2NoSpike_depthNorm.bw ${workDir}/bigwig/H3K27ac_EpiLC_Rep1_depthNorm.bw ${workDir}/bigwig/H3K27ac_EpiLC_Rep2_depthNorm.bw \
--beforeRegionStartLength 2000 \
--afterRegionStartLength 2000 \
--binSize 25 \
--missingDataAsZero \
-p 8 \
-o ${outDir}/matrix_diff_H3K27ac_onlyChange.gz

plotHeatmap \
-m ${outDir}/matrix_diff_H3K27ac_onlyChange.gz \
--heatmapHeight 0.5 \
--heatmapWidth 3 \
--colorMap Blues Blues Greens Greens \
-T "" --refPointLabel "" \
--samplesLabel "ESC1" "ESC2" "EpiLC1" "EpiLC2" \
--legendLocation none \
--yAxisLabel "RPM" \
--xAxisLabel "peak center" \
--regionsLabel "ESC" "EpiLC" \
-o ${outDir}/heatmap_diff_H3K27ac_onlyChange.svg

# REVIEWED FOR PUBLICATION END
############################################################################################



# Extracting RAD21 specific peaks
mamba activate bioinfo
chipDir="/Volumes/UKJIN_SSD/reference/ChIP/RAD21-HA"
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


BW="${chipDir}/RAD21_ESC-DMSO_Rep1_spikeNorm_inputNorm.bw"
BASENAME=$(basename "${BED1%.*}")

# 2) Convert bigWig to bedGraph
bigWigToBedGraph "$BW" "$OUTDIR/${BASENAME}.signal.bedgraph"

# 3) Map bigWig signal over peaks:
#    Append two columns to the BED: max and mean signal per peak
bedtools map \
  -a "$BED1" \
  -b "$OUTDIR/${BASENAME}.signal.bedgraph" \
  -c 4 -o max,mean \
  > "$OUTDIR/${BASENAME}.DMSO.with_signal.bed"

# 4) Sort by max signal descending and take top N
#    The two appended columns are the last two fields; we sort on the penultimate field (max)
awk -v OFS='\t' '{print $0, (NF-1), NF}' "$OUTDIR/${BASENAME}.DMSO.with_signal.bed" | \
  sort -k$(awk 'END{print 0}' /dev/null) >/dev/null 2>&1 # placeholder to avoid shellcheck
# robust sort without awk gymnastics:
# extract max column index, then sort
MAXCOL=$(head -n1 "$OUTDIR/${BASENAME}.DMSO.with_signal.bed" | awk '{print NF-1}')
sort -k"${MAXCOL}","${MAXCOL}"gr \
  "$OUTDIR/${BASENAME}.DMSO.with_signal.bed" \
  | head -n "$TOPN" \
  > "$OUTDIR/${BASENAME}.top${TOPN}_by_max_signal.bed"
