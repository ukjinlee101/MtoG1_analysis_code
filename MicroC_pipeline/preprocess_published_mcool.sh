conda activate hicexplorer
COOLDIR="data"


EXP="GSE178982_AsyncUT_pooled"
EXP="GSE178982_AsyncAID_pooled"

# GSE178982 has weight column calculated using ICE normalization
for RES in 1000 5000 10000 100000; do
    IN_RES="${COOLDIR}/${EXP}_allRes.mcool::/resolutions/${RES}"
    KR_OUT="${COOLDIR}/${EXP}_${RES}bp_KR.cool"

    # Extract raw .cool at each resolution
    hicConvertFormat \
        -m "${IN_RES}" \
        --inputFormat cool \
        --outputFormat cool \
        -o "${KR_OUT}"
done

# Make .cool in 25kb resolution by coarsening the 5kb
COARSEN_IN="${COOLDIR}/${EXP}_5000bp_KR.cool"
COARSEN_TMP="${COOLDIR}/${EXP}_25000bp_KR.cool"
cooler coarsen -k 5 -p 32 -o "${COARSEN_TMP}" "${COARSEN_IN}"

# Calculate KR normalization and expected reads.
for RES in 1000 5000 10000 25000 100000; do
    KR_OUT="${COOLDIR}/${EXP}_${RES}bp_KR.cool"
    EXP_OUT="${COOLDIR}/${EXP}_${RES}bp_KR_exp.tsv"

    echo "Balancing (KR) at ${RES}bp..."
    # recompute KR weights into bins.weight
    cooler balance --force -p 16 \
        "${KR_OUT}"

    echo "Computing expected-cis at ${RES}bp..."
    # compute expected‐cis
    cooltools expected-cis \
        -p 16 \
        -o "${EXP_OUT}" \
        "${KR_OUT}"
done