conda activate chromosight

# Quantify chromosight scores based on Hansen loops

sample=("GSE178982_AsyncUT_pooled" "GSE178982_AsyncAID_pooled" "G1DMSO_pooled" "G1dTAG_pooled" "G1A485_pooled" "EpiG1DMSO_pooled" "EpiG1dTAG_pooled")


resolutions=("20000" "10000" "4000" "2000" "1000" "800" "600" "400")


coolDir="/athena/apostoloulab/scratch/ukl4001/data/cool_norm_pooled"
resultDir="/athena/apostoloulab/scratch/ukl4001/data/loop_Tjian"
mkdir -p $resultDir

for res in "${resolutions[@]}"; do
    for samp in "${sample[@]}"; do
        chromosight quantify ${resultDir}/TjianChromo_union_${res}bp.bedpe \
                        ${coolDir}/${samp}_${res}bp_KR.cool \
                        ${resultDir}/${samp}_union_${res}bp_pu100pz100 \
                        --pattern=loops \
                        --threads=8 \
                        --perc-undetected=100 \
                        --perc-zero=100
    done
done
