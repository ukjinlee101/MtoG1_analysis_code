mamba activate bioinfo

### Paired-end reads is used to make virtual 4C
### Converting validPairs (.pairs) into BEDPE (chr/start/end/chr/start/end/*/*/+/-) format
##############################################################################
PAIRDIR="/athena/apostoloulab/scratch/ukl4001/temp"
REFDIR="/athena/apostoloulab/scratch/ukl4001/reference"
GENOMEDIR="/athena/apostoloulab/scratch/collab/Genome/Mus_Musculus/UCSC/mm10/Sequence/WholeGenomeFasta"
VPDIR="/athena/apostoloulab/scratch/ukl4001/MtoG1_analysis_code/MicroC_v4C/vpInfo"
TEMPDIR="/athena/apostoloulab/scratch/ukl4001/temp"
AIDIR="/athena/apostoloulab/scratch/ukl4001/data/v4C/allInfo"
mkdir -p ${AIDIR}

### Splitting per chr
mkdir -p ${PAIRDIR}/mappedPairChr
cd ${PAIRDIR}/mappedPairChr

for EXP in "G1DMSO_pooled" "G1dTAG_pooled" "G1A485_pooled" "EpiG1DMSO_pooled" "EpiG1dTAG_pooled" "GSE178982_AsyncUT_pooled" "GSE178982_AsyncAID_pooled"
do
    grep -v "^#" ${PAIRDIR}/mappedPair/${EXP}.pairs > ${PAIRDIR}/temp.${EXP}
    awk -F'\t' '{print > ("'${PAIRDIR}/mappedPairChr/${EXP}'."$2)}' ${PAIRDIR}/temp.${EXP}
    rm ${PAIRDIR}/temp.${EXP}


    for chr in "chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chrX" "chrY" "chrM"
    do
        awk '{print $2"\t"$3"\t"$3+100"\t"$4"\t"$5"\t"$5+100}' ${PAIRDIR}/mappedPairChr/${EXP}.${chr} > ${PAIRDIR}/mappedPairChr/${EXP}.${chr}.bedpe
        rm ${PAIRDIR}/mappedPairChr/${EXP}.${chr}
    done

    touch ${PAIRDIR}/${EXP}.DONE
done


#######################################################################################


# Make vp info file per gene from vp.all.bed file
awk -F'\t' '{print $2"\t"$3"\t"$3+1 > ("'${VPDIR}'/vp."$1".bed")}' ${VPDIR}/vp.all.bed
##############################################################################
# MAKE REGION OF INTEREST FOR EACH GENE 10kb
for gene in "Fosl2" "Myc" "Klf4" "Tbx3"
do
    # Extract bin with VP (region1)
    intersectBed -a ${REFDIR}/mm10/mm10.bin.10kb.bed -b ${VPDIR}/vp.${gene}.bed -wa -u > ${TEMPDIR}/region1.${gene}
    # Expand the region of interest to +- 2Mb (slop)
    slopBed -i ${TEMPDIR}/region1.${gene} -g ${GENOMEDIR}/genome.fa.fai -b 2000000 > ${TEMPDIR}/slop.${gene}
    # Create 10kb windows with 500 bp steps within 4Mb regions (step)
    bedtools makewindows -b ${TEMPDIR}/slop.${gene} -w 10000 -s 500 > ${TEMPDIR}/step.${gene}
    # Count number of restriction sites within each window (step1)
    intersectBed -a ${TEMPDIR}/step.${gene} -b ${REFDIR}/mm10/MboI_mm10.bed -c > ${TEMPDIR}/step1.${gene}

    rm ${TEMPDIR}/slop.${gene} ${TEMPDIR}/step.${gene}
done



for EXP in "G1DMSO_pooled" "G1dTAG_pooled" "G1A485_pooled" "EpiG1DMSO_pooled" "EpiG1dTAG_pooled" "GSE178982_AsyncUT_pooled" "GSE178982_AsyncAID_pooled"
do
# Extract read pairs which has one anchor at the 10 kb bin where viewpoint is located
    for gene in "Fosl2" "Myc" "Klf4" "Tbx3"
    do
        echo "Processing EXP: ${EXP}, Gene: ${gene}"
        firstline=$(head -n 1 ${VPDIR}/vp.${gene}.bed)
        chr=$(echo ${firstline} | awk '{print $1}')

        pairToBed -a ${PAIRDIR}/mappedPairChr/${EXP}.${chr}.bedpe -b ${TEMPDIR}/region1.${gene} -type "either" > ${TEMPDIR}/pairs.${gene}
        cut -f1,2,3 ${TEMPDIR}/pairs.${gene} > ${TEMPDIR}/left.${gene}
        cut -f4,5,6 ${TEMPDIR}/pairs.${gene} > ${TEMPDIR}/right.${gene}
        cat ${TEMPDIR}/left.${gene} ${TEMPDIR}/right.${gene} > ${TEMPDIR}/singlereads.${gene}
        #calculate read coverage for all bins around the region of interest
        intersectBed -a ${TEMPDIR}/step1.${gene} -b ${TEMPDIR}/singlereads.${gene} -c > ${TEMPDIR}/enrichment.${gene}
        # add location of region of interest - We will delete this from downstream analysis for plots - it is useful to keep it as a separate column.
        intersectBed -a ${TEMPDIR}/enrichment.${gene} -b ${TEMPDIR}/region1.${gene} -c > ${AIDIR}/allinfo.10kb.${EXP}.${gene}
        rm ${TEMPDIR}/pairs.${gene} ${TEMPDIR}/left.${gene} ${TEMPDIR}/right.${gene} ${TEMPDIR}/singlereads.${gene}  ${TEMPDIR}/enrichment.${gene}
    done
done

##############################################################################
# MAKE REGION OF INTEREST FOR EACH GENE 5kb
for gene in "Fosl2" "Myc" "Klf4" "Tbx3"
do
    # Extract bin with VP (region1)
    #intersectBed -a ${REFDIR}/mm10/mm10.bin.5kb.bed -b ${VPDIR}/vp.${gene}.bed -wa -u > ${TEMPDIR}/region1.${gene}
    # Expand the region of interest to +- 2Mb (slop)
    slopBed -i ${TEMPDIR}/region1.${gene} -g ${GENOMEDIR}/genome.fa.fai -b 2000000 > ${TEMPDIR}/slop.${gene}
    # Create 5kb windows with 500 bp steps within 4Mb regions (step)
    bedtools makewindows -b ${TEMPDIR}/slop.${gene} -w 5000 -s 500 > ${TEMPDIR}/step.${gene}
    # Count number of restriction sites within each window (step1)
    intersectBed -a ${TEMPDIR}/step.${gene} -b ${REFDIR}/mm10/MboI_mm10.bed -c > ${TEMPDIR}/step1.${gene}

    rm ${TEMPDIR}/slop.${gene} ${TEMPDIR}/step.${gene}
done



for EXP in "G1DMSO_pooled" "G1dTAG_pooled" "G1A485_pooled" "EpiG1DMSO_pooled" "EpiG1dTAG_pooled" "GSE178982_AsyncUT_pooled" "GSE178982_AsyncAID_pooled"
do
# Extract read pairs which has one anchor at the 10 kb bin where viewpoint is located
    for gene in "Fosl2" "Myc" "Klf4" "Tbx3"
    do
        echo "Processing EXP: ${EXP}, Gene: ${gene}"
        firstline=$(head -n 1 ${VPDIR}/vp.${gene}.bed)
        chr=$(echo ${firstline} | awk '{print $1}')

        pairToBed -a ${PAIRDIR}/mappedPairChr/${EXP}.${chr}.bedpe -b ${TEMPDIR}/region1.${gene} -type "either" > ${TEMPDIR}/pairs.${gene}
        cut -f1,2,3 ${TEMPDIR}/pairs.${gene} > ${TEMPDIR}/left.${gene}
        cut -f4,5,6 ${TEMPDIR}/pairs.${gene} > ${TEMPDIR}/right.${gene}
        cat ${TEMPDIR}/left.${gene} ${TEMPDIR}/right.${gene} > ${TEMPDIR}/singlereads.${gene}
        #calculate read coverage for all bins around the region of interest
        intersectBed -a ${TEMPDIR}/step1.${gene} -b ${TEMPDIR}/singlereads.${gene} -c > ${TEMPDIR}/enrichment.${gene}
        # add location of region of interest - We will delete this from downstream analysis for plots - it is useful to keep it as a separate column.
        intersectBed -a ${TEMPDIR}/enrichment.${gene} -b ${TEMPDIR}/region1.${gene} -c > ${AIDIR}/allinfo.5kb.${EXP}.${gene}
        rm ${TEMPDIR}/pairs.${gene} ${TEMPDIR}/left.${gene} ${TEMPDIR}/right.${gene} ${TEMPDIR}/singlereads.${gene}  ${TEMPDIR}/enrichment.${gene}
    done
done


##############################################################################
# Checking stats of pairs file for normalization
mamba activate MicroC_QC
STATDIR="/athena/apostoloulab/scratch/ukl4001/data/v4C/stats"
mkdir -p ${STATDIR}
for EXP in "G1DMSO_pooled" "G1dTAG_pooled" "G1A485_pooled" "EpiG1DMSO_pooled" "EpiG1dTAG_pooled" "GSE178982_AsyncUT_pooled" "GSE178982_AsyncAID_pooled"
do
    echo "Running pairtools stats on ${EXP}"
    pairtools stats ${PAIRDIR}/mappedPair/${EXP}.pairs > ${STATDIR}/${EXP}.stats.txt
done
