# DRAWING METAPLOT & TORNADO PLOT ON UNION RAD21 PEAKS TO SHOW DEPLETION UPON dTAG TREATMENT
mamba activate macs2
chipDir="/Volumes/UKJIN_SSD/reference/ChIP"

# Merge peaks from ESC DMSO and ESC dTAG
cat ${chipDir}/RAD21_ESC-DMSO_unionPeaks.bed ${chipDir}/RAD21_ESC-dTAG_unionPeaks.bed \
| bedtools sort -i - \
| bedtools merge -d 200 -i - \
> "${chipDir}/RAD21_ESC_DMSOdTAG_unionPeaks.bed"


dataDir="/Volumes/UKJIN_SSD/data/deeptools"
figDir=""/Volumes/UKJIN_SSD/figure/deeptools""
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
--heatmapHeight 6 \
--colorMap Blues -max 1 \
-T $name --refPointLabel "" \
--samplesLabel "DMSO1" "DMSO2" "dTAG1" "dTAG2" \
--legendLocation none \
--yAxisLabel "spikeInputNorm" \
--xAxisLabel "peak center" \
--plotFileFormat svg \
-o ${figDir}/heatmap_${name}.svg
