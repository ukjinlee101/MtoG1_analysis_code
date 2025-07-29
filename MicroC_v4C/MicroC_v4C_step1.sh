	conda activate bioinfo


### Paired-end reads is used to make virtual 4C
### Converting validPairs (.pairs) into BEDPE (chr/start/end/chr/start/end/*/*/+/-) format
##############################################################################
REFDIR="/Volumes/UKJIN_SSD/Genomics_03_Analysis_Working/data_vault_2024summer_microC/reference"
GENOMEDIR="/Volumes/Charmander/Genomics_00_Reference/Genome/M.musculus/mm10/WholeGenomeFasta"
VPDIR="/Volumes/UKJIN_SSD/Genomics_03_Analysis_Working/data_vault_2024summer_microC/source/v4C/vpInfo"
TEMPDIR="/Volumes/UKJIN_SSD/Genomics_03_Analysis_Working/data_vault_2024summer_microC/source/v4C/temp"
AIDIR="/Volumes/UKJIN_SSD/Genomics_03_Analysis_Working/data_vault_2024summer_microC/source/v4C/allinfo"
PAIRDIR="/Volumes/Charmander/PASSED_Genomics_02_Pipeline_Archive/temp"

mkdir ${PAIRDIR}/mappedPairChr
mkdir ${TEMPDIR}
mkdir ${AIDIR}
##############################################################################
### Splitting per chr
cd ${PAIRDIR}/mappedPairChr
EXP="G1.dTAG.Merged"
EXP="G1.A485.Merged"
EXP="G1.DMSO.Merged"

for EXP in "G1.A485.Merged" "G1.DMSO.Merged"
do
	grep -v "^#" ${PAIRDIR}/mappedPair/${EXP}.pairs > ${PAIRDIR}/temp.${EXP}
	awk -F'\t' '{print > ("'${PAIRDIR}/mappedPairChr/${EXP}'."$2)}' ${PAIRDIR}/temp.${EXP}
	rm ${PAIRDIR}/temp.${EXP}


	for chr in "chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chrX" "chrY" "chrM"
	do
		awk '{print $2"\t"$3"\t"$3+100"\t"$4"\t"$5"\t"$5+100}' ${PAIRDIR}/mappedPairChr/${EXP}.${chr} > ${PAIRDIR}/mappedPairChr/${EXP}.${chr}.bedpe
		rm ${PAIRDIR}/mappedPairChr/${EXP}.${chr}
	done
done
##############################################################################
# Make vp info file per gene from vp.all.bed file
awk -F'\t' '{print $2"\t"$3"\t"$3+1 > ("'${VPDIR}'/vp."$1".bed")}' ${VPDIR}/vp.all.bed
##############################################################################
# MAKE REGION OF INTEREST FOR EACH GENE
for gene in "Cited2" "Dusp4" "Dusp5" "Epha2" "Evx1" "Fosl1" "Fosl2" "Jun" "Klf4" "Myc" "Phlda1" "Sall1" "Spry4" "Tbx3" "Wnt3" "Wnt6" "Xbp1"
do
	# Extract bin with VP (region1)
	intersectBed -a ${REFDIR}/mm10.bin.10kb.bed -b ${VPDIR}/vp.${gene}.bed -wa -u > ${TEMPDIR}/region1.${gene}
	# Expand the region of interest to +- 2Mb (slop)
	slopBed -i ${TEMPDIR}/region1.${gene} -g ${GENOMEDIR}/genome.fa.fai -b 2000000 > ${TEMPDIR}/slop.${gene}
	# Create 10kb windows with 500 bp steps within 4Mb regions (step)
	bedtools makewindows -b ${TEMPDIR}/slop.${gene} -w 10000 -s 500 > ${TEMPDIR}/step.${gene}
	# Count number of restriction sites within each window (step1)
	intersectBed -a ${TEMPDIR}/step.${gene} -b ${REFDIR}/MboI_mm10.bed -c > ${TEMPDIR}/step1.${gene}

	rm ${TEMPDIR}/slop.${gene} ${TEMPDIR}/step.${gene}
done

# 2024.08.26
for gene in "Trib1" "Olig2" "Fam43a" "Epha4" "1700017B05Rik" "Gabbr1" "Rnd3" "Irs2" "Cxadr" "Hs6st1"
do
	# Extract bin with VP (region1)
	intersectBed -a ${REFDIR}/mm10.bin.10kb.bed -b ${VPDIR}/vp.${gene}.bed -wa -u > ${TEMPDIR}/region1.${gene}
	# Expand the region of interest to +- 2Mb (slop)
	slopBed -i ${TEMPDIR}/region1.${gene} -g ${GENOMEDIR}/genome.fa.fai -b 2000000 > ${TEMPDIR}/slop.${gene}
	# Create 10kb windows with 500 bp steps within 4Mb regions (step)
	bedtools makewindows -b ${TEMPDIR}/slop.${gene} -w 10000 -s 500 > ${TEMPDIR}/step.${gene}
	# Count number of restriction sites within each window (step1)
	intersectBed -a ${TEMPDIR}/step.${gene} -b ${REFDIR}/MboI_mm10.bed -c > ${TEMPDIR}/step1.${gene}

	rm ${TEMPDIR}/slop.${gene} ${TEMPDIR}/step.${gene}
done
##############################################################################
 
for EXP in "G1.DMSO.Merged" "G1.A485.Merged" "G1.dTAG.Merged"
do
# Extract read pairs which has one anchor at the 10 kb bin where viewpoint is located
	for gene in "Cited2" "Dusp4" "Dusp5" "Epha2" "Evx1" "Fosl1" "Fosl2" "Jun" "Klf4" "Myc" "Phlda1" "Sall1" "Spry4" "Tbx3" "Wnt3" "Wnt6" "Xbp1"
	do
	firstline=$(head -n 1 ${VPDIR}/vp.${gene}.bed)
	chr=$(echo ${firstline} | awk '{print $1}')

	pairToBed -a ${PAIRDIR}/mappedPairChr/${EXP}.${chr}.bedpe -b ${TEMPDIR}/region1.${gene} -type "either" > ${TEMPDIR}/pairs.${gene}
	cut -f1,2,3 ${TEMPDIR}/pairs.${gene} > ${TEMPDIR}/left.${gene}
	cut -f4,5,6 ${TEMPDIR}/pairs.${gene} > ${TEMPDIR}/right.${gene}
	cat ${TEMPDIR}/left.${gene} ${TEMPDIR}/right.${gene} > ${TEMPDIR}/singlereads.${gene}
	#calculate read coverage for all bins around the region of interest
	intersectBed -a ${TEMPDIR}/step1.${gene} -b ${TEMPDIR}/singlereads.${gene} -c > ${TEMPDIR}/enrichment.${gene}
	# add location of region of interest - We will delete this from downstream analysis for plots - it is useful to keep it as a separate column.
	intersectBed -a ${TEMPDIR}/enrichment.${gene} -b ${TEMPDIR}/region1.${gene} -c > ${AIDIR}/allinfo.${EXP}.${gene}
	rm ${TEMPDIR}/pairs.${gene} ${TEMPDIR}/left.${gene} ${TEMPDIR}/right.${gene} ${TEMPDIR}/singlereads.${gene}  ${TEMPDIR}/enrichment.${gene}
	done
done

# 2024.08.26

for EXP in "G1.DMSO.Merged" "G1.A485.Merged" "G1.dTAG.Merged"
do
# Extract read pairs which has one anchor at the 10 kb bin where viewpoint is located
	for gene in "Trib1" "Olig2" "Fam43a" "Epha4" "1700017B05Rik" "Gabbr1" "Rnd3" "Irs2" "Cxadr" "Hs6st1"
	do
	firstline=$(head -n 1 ${VPDIR}/vp.${gene}.bed)
	chr=$(echo ${firstline} | awk '{print $1}')

	pairToBed -a ${PAIRDIR}/mappedPairChr/${EXP}.${chr}.bedpe -b ${TEMPDIR}/region1.${gene} -type "either" > ${TEMPDIR}/pairs.${gene}
	cut -f1,2,3 ${TEMPDIR}/pairs.${gene} > ${TEMPDIR}/left.${gene}
	cut -f4,5,6 ${TEMPDIR}/pairs.${gene} > ${TEMPDIR}/right.${gene}
	cat ${TEMPDIR}/left.${gene} ${TEMPDIR}/right.${gene} > ${TEMPDIR}/singlereads.${gene}
	#calculate read coverage for all bins around the region of interest
	intersectBed -a ${TEMPDIR}/step1.${gene} -b ${TEMPDIR}/singlereads.${gene} -c > ${TEMPDIR}/enrichment.${gene}
	# add location of region of interest - We will delete this from downstream analysis for plots - it is useful to keep it as a separate column.
	intersectBed -a ${TEMPDIR}/enrichment.${gene} -b ${TEMPDIR}/region1.${gene} -c > ${AIDIR}/allinfo.${EXP}.${gene}
	rm ${TEMPDIR}/pairs.${gene} ${TEMPDIR}/left.${gene} ${TEMPDIR}/right.${gene} ${TEMPDIR}/singlereads.${gene}  ${TEMPDIR}/enrichment.${gene}
	done
done

##############################################################################
# Checking stats of pairs file for normalization
for EXP in "G1.DMSO.Merged" "G1.A485.Merged" "G1.dTAG.Merged"
do
	pairtools stats -o ${PAIRDIR}/${EXP}.stats.txt ${PAIRDIR}/mappedPair/${EXP}.pairs
done

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
# REPEATING FOR 5KB
# MAKE REGION OF INTEREST FOR EACH GENE
for gene in "Cited2" "Dusp4" "Dusp5" "Epha2" "Evx1" "Fosl1" "Fosl2" "Jun" "Klf4" "Myc" "Phlda1" "Sall1" "Spry4" "Tbx3" "Wnt3" "Wnt6" "Xbp1"
do
	# Extract bin with VP (region1)
	intersectBed -a ${REFDIR}/mm10.bin.5kb.bed -b ${VPDIR}/vp.${gene}.bed -wa -u > ${TEMPDIR}/region1.${gene}
	# Expand the region of interest to +- 2Mb (slop)
	slopBed -i ${TEMPDIR}/region1.${gene} -g ${GENOMEDIR}/genome.fa.fai -b 2000000 > ${TEMPDIR}/slop.${gene}
	# Create 10kb windows with 500 bp steps within 4Mb regions (step)
	bedtools makewindows -b ${TEMPDIR}/slop.${gene} -w 5000 -s 500 > ${TEMPDIR}/step.${gene}
	# Count number of restriction sites within each window (step1)
	intersectBed -a ${TEMPDIR}/step.${gene} -b ${REFDIR}/MboI_mm10.bed -c > ${TEMPDIR}/step1.${gene}

	rm ${TEMPDIR}/slop.${gene} ${TEMPDIR}/step.${gene}
done

# 2024.08.26
for gene in "Trib1" "Olig2" "Fam43a" "Epha4" "1700017B05Rik" "Gabbr1" "Rnd3" "Irs2" "Cxadr" "Hs6st1"
do
	# Extract bin with VP (region1)
	intersectBed -a ${REFDIR}/mm10.bin.5kb.bed -b ${VPDIR}/vp.${gene}.bed -wa -u > ${TEMPDIR}/region1.${gene}
	# Expand the region of interest to +- 2Mb (slop)
	slopBed -i ${TEMPDIR}/region1.${gene} -g ${GENOMEDIR}/genome.fa.fai -b 2000000 > ${TEMPDIR}/slop.${gene}
	# Create 10kb windows with 500 bp steps within 4Mb regions (step)
	bedtools makewindows -b ${TEMPDIR}/slop.${gene} -w 5000 -s 500 > ${TEMPDIR}/step.${gene}
	# Count number of restriction sites within each window (step1)
	intersectBed -a ${TEMPDIR}/step.${gene} -b ${REFDIR}/MboI_mm10.bed -c > ${TEMPDIR}/step1.${gene}

	rm ${TEMPDIR}/slop.${gene} ${TEMPDIR}/step.${gene}
done
##############################################################################
 
for EXP in "G1.DMSO.Merged" "G1.A485.Merged" "G1.dTAG.Merged"
do
# Extract read pairs which has one anchor at the 5 kb bin where viewpoint is located
	for gene in "Cited2" "Dusp4" "Dusp5" "Epha2" "Evx1" "Fosl1" "Fosl2" "Jun" "Klf4" "Myc" "Phlda1" "Sall1" "Spry4" "Tbx3" "Wnt3" "Wnt6" "Xbp1"
	do
	firstline=$(head -n 1 ${VPDIR}/vp.${gene}.bed)
	chr=$(echo ${firstline} | awk '{print $1}')

	pairToBed -a ${PAIRDIR}/mappedPairChr/${EXP}.${chr}.bedpe -b ${TEMPDIR}/region1.${gene} -type "either" > ${TEMPDIR}/pairs.${gene}
	cut -f1,2,3 ${TEMPDIR}/pairs.${gene} > ${TEMPDIR}/left.${gene}
	cut -f4,5,6 ${TEMPDIR}/pairs.${gene} > ${TEMPDIR}/right.${gene}
	cat ${TEMPDIR}/left.${gene} ${TEMPDIR}/right.${gene} > ${TEMPDIR}/singlereads.${gene}
	#calculate read coverage for all bins around the region of interest
	intersectBed -a ${TEMPDIR}/step1.${gene} -b ${TEMPDIR}/singlereads.${gene} -c > ${TEMPDIR}/enrichment.${gene}
	# add location of region of interest - We will delete this from downstream analysis for plots - it is useful to keep it as a separate column.
	intersectBed -a ${TEMPDIR}/enrichment.${gene} -b ${TEMPDIR}/region1.${gene} -c > ${AIDIR}/allinfo-5kb.${EXP}.${gene}
	rm ${TEMPDIR}/pairs.${gene} ${TEMPDIR}/left.${gene} ${TEMPDIR}/right.${gene} ${TEMPDIR}/singlereads.${gene}  ${TEMPDIR}/enrichment.${gene}
	done
done

# 2024.08.26
for EXP in "G1.DMSO.Merged" "G1.A485.Merged" "G1.dTAG.Merged"
do
# Extract read pairs which has one anchor at the 5 kb bin where viewpoint is located
	for gene in "Trib1" "Olig2" "Fam43a" "Epha4" "1700017B05Rik" "Gabbr1" "Rnd3" "Irs2" "Cxadr" "Hs6st1"
	do
	firstline=$(head -n 1 ${VPDIR}/vp.${gene}.bed)
	chr=$(echo ${firstline} | awk '{print $1}')

	pairToBed -a ${PAIRDIR}/mappedPairChr/${EXP}.${chr}.bedpe -b ${TEMPDIR}/region1.${gene} -type "either" > ${TEMPDIR}/pairs.${gene}
	cut -f1,2,3 ${TEMPDIR}/pairs.${gene} > ${TEMPDIR}/left.${gene}
	cut -f4,5,6 ${TEMPDIR}/pairs.${gene} > ${TEMPDIR}/right.${gene}
	cat ${TEMPDIR}/left.${gene} ${TEMPDIR}/right.${gene} > ${TEMPDIR}/singlereads.${gene}
	#calculate read coverage for all bins around the region of interest
	intersectBed -a ${TEMPDIR}/step1.${gene} -b ${TEMPDIR}/singlereads.${gene} -c > ${TEMPDIR}/enrichment.${gene}
	# add location of region of interest - We will delete this from downstream analysis for plots - it is useful to keep it as a separate column.
	intersectBed -a ${TEMPDIR}/enrichment.${gene} -b ${TEMPDIR}/region1.${gene} -c > ${AIDIR}/allinfo-5kb.${EXP}.${gene}
	rm ${TEMPDIR}/pairs.${gene} ${TEMPDIR}/left.${gene} ${TEMPDIR}/right.${gene} ${TEMPDIR}/singlereads.${gene}  ${TEMPDIR}/enrichment.${gene}
	done
done

