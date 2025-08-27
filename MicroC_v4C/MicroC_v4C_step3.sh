mamba activate bioinfo

V4CDIR="/Volumes/UKJIN_SSD/data/v4C"
REFDIR="/Volumes/UKJIN_SSD/reference/mm10"

mkdir ${V4CDIR}/bw
##############################################################################
for file in ${V4CDIR}/bedgraph/*.bedGraph;
do
	filename=$(basename "$file" .bedGraph)
	echo $filename

	sort -k1,1 -k2,2n $file > ${V4CDIR}/bw/temp
	bedGraphToBigWig ${V4CDIR}/bw/temp ${REFDIR}/mm10.chrom.sizes ${V4CDIR}/bw/${filename}.bw
	rm ${V4CDIR}/bw/temp
done

