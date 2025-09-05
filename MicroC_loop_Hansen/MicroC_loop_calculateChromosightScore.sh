conda activate chromosight

# Quantify chromosight scores based on Hansen loops

sample=("GSE178982_AsyncUT_pooled" "GSE178982_AsyncAID_pooled" "G1DMSO_pooled" "G1dTAG_pooled" "G1A485_pooled" "EpiG1DMSO_pooled" "EpiG1dTAG_pooled")
resolutions=("5000" "2000" "1000")

coolDir="/athena/apostoloulab/scratch/ukl4001/data/cool_norm_pooled"
resultDir="/athena/apostoloulab/scratch/ukl4001/data/loop_Hansen"
mkdir -p $resultDir

for res in "${resolutions[@]}"; do
    for samp in "${sample[@]}"; do
        chromosight quantify ${resultDir}/Hansen_union_allRes_postprocessed_${res}bp.bedpe \
                        ${coolDir}/${samp}_${res}bp_KR.cool \
                        ${resultDir}/${samp}_union_${res}bp_pu100pz100 \
                        --pattern=loops \
                        --threads=32 \
                        --perc-undetected=100 \
                        --perc-zero=100
    done
done
