# CONVERT VALIDPAIRS TO hicSummary

PAIRDIR="/athena/apostoloulab/scratch/ukl4001/MicroC_MtoG1_Batch123_pooled/pipeline/result/mappedPair_pooled"
PAIRDIR="/athena/apostoloulab/scratch/ukl4001/MicroC_Async_GSE178982/pipeline/result/mappedPair_pooled"
OUTDIR="/athena/apostoloulab/scratch/ukl4001/temp"

sample="G1DMSO_pooled"
sample="G1dTAG_pooled"
sample="G1A485_pooled"
sample="GSE178982_AsyncUT_pooled"
sample="GSE178982_AsyncAID_pooled"


infile="${PAIRDIR}/${sample}.pairs"
tempfile="${OUTDIR}/${sample}.output.txt"
outfile="${OUTDIR}/hicSummary/${sample}.summary"

# Strip header up to the line matching "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type"
sed '1,/#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type/d' \
    "${infile}" > "${tempfile}"

# Reorder columns: readID, chrom1, pos1, strand1, chrom2, pos2, pair_type
awk -F $'\t' 'BEGIN {OFS = FS} { print $1, $2, $3, $6, $4, $5, $7 }' "${tempfile}" > "${outfile}"

rm ${tempfile}

# MAKE 10kb BIN
SCRDIR="../script"
REFDIR="../reference/mm10"
"$SCRDIR"/binGenome.py  ${REFDIR}/mm10.chrom.sizes.canonical 10000 ${REFDIR}/mm10.bin.10kb.bed.gz
gzip -d ${REFDIR}/mm10.bin.10kb.bed.gz


# RUN CscoreTool1.1
# REFERENCE
# https://github.com/scoutzxb/CscoreTool

# RUN
# a. windows.bed This file is to specify the genomic windows to analyze. It should be equal-length bed files covering the region of interest, presumably a whole chromosome or whole genome. An example hg19_chr1_10k.bed can be downloaded.These files can also be generatd using the generateEqualLengthBed.cpp program provided here. The chromosome size files used for generating windows.bed are also available for download.

# b. input.summary This file is the main input file for Hi-C interactions. We accept the same format as the HiCsummary file format for HOMER runHiCpca.pl. See http://homer.ucsd.edu/homer/interactions/HiCtagDirectory.html An example file test.summary.gz can be downloaded. This is 0.5% randomly selected reads in chr1 from the High-resolution GM12878 cell Hi-C dataset (Rao, 2014).

# c. OutputPrefix This is the prefix for output files.

# d. session This the number of sessions to use. The number of choide depends on the resource available. The program is not perfectly parallelized, so using more sessions just partly improve the speed.

# e. minDis This is the criteria of minimum interaction distance to be considered. It should be no less than the analysis resolution. We suggest using 1000000 (1M), because interactions of shorter distance may be more affected by TAD structure than A/B compartments.

# f. chrName This is a parameter by choice. Users can speficy a certain chromosome to analyze using it.

sbatch cscore_G1DMSO.sh
sbatch cscore_G1dTAG.sh
sbatch cscore_G1A485.sh

sbatch cscore_GSE178982_AsyncUT.sh
sbatch cscore_GSE178982_AsyncAID.sh