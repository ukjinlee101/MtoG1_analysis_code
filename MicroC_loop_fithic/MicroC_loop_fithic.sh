# FitHiC2 requires 2 inputs. Contact founts and fragments
mamba activate fithic2
cd MtoG1_analysis_code/MicroC_loop_fithic/

# Making fragments
./createFitHiCFragments-fixedsize.py \
    --chrLens ../../reference/mm10/mm10.chrom.sizes \
    --outFile ../../data/loop_fithic/fragments.mm10.10kb.txt.gz \
    --resolution 10000


# Making interactions from hic & calculate bias for KR normalization
for sample in G1DMSO_pooled G1dTAG_pooled G1A485_pooled EpiG1DMSO_pooled EpiG1dTAG_pooled; do
    echo "Processing sample: $sample"
    DATADIR="../../data/loop_fithic"
    OUTDIR="${DATADIR}/${sample}"
    mkdir -p ${OUTDIR}

    MERGED_INTERACTIONS_FILE="${OUTDIR}/interactions_${sample}_10kb.txt"
    > ${MERGED_INTERACTIONS_FILE}

    for chr in chr{1..19} chrX; do
        TEMP_OUT_FILE="../../data/loop_fithic/interactions_${sample}_${chr}_${chr}_10kb.txt"

        ./createFitHiCContacts-hic.py \
            --HiCFile ../../data/hic_pooled/${sample}_allRes.hic \
            --CHR1 ${chr} \
            --CHR2 ${chr} \
            --resolution 10000 \
            --Norm "NONE" \
            --datatype observed \
            --outFile ${TEMP_OUT_FILE}

        cat ${TEMP_OUT_FILE} >> ${MERGED_INTERACTIONS_FILE}
        rm ${TEMP_OUT_FILE}
    done

    echo "Gzipping merged file for sample: $sample"    
    gzip ${MERGED_INTERACTIONS_FILE}

    ./HiCKRy.py \
        -i ${MERGED_INTERACTIONS_FILE}.gz \
        -f ${DATADIR}/fragments.mm10.10kb.txt.gz \
        -o ${OUTDIR}/bias_${sample}_10kb_0.1.txt.gz \
        -x 0.1
done


# Running FitHiC2
sample="EpiG1DMSO_pooled"
DATADIR="../../data/loop_fithic"
OUTDIR="${DATADIR}/${sample}"
mkdir -p ${OUTDIR}
fithic -i ${OUTDIR}/interactions_${sample}_10kb.txt.gz \
    -f ${DATADIR}/fragments.mm10.10kb.txt.gz \
    -o ${OUTDIR}/ \
    -r 10000 \
    -t ${OUTDIR}/bias_${sample}_10kb_0.1.txt.gz \
    -v