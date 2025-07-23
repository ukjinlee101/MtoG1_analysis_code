#!/bin/bash

genome=mm10
destTmp="./tmp"
destRegs="./regions"
mkdir -p $destTmp
mkdir -p $destRegs
#wget --directory-prefix=$dest http://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz


if [ ! -f  $destTmp/ucsc.tables.mod.txt ] ; then	
	echo "show tables;" | mysql --user=genome --host=genome-mysql.cse.ucsc.edu -D $genome -A > $destTmp/ucsc.tables.txt
	
	echo "please select tables in $destTmp/ucsc.tables.txt and save as $destTmp/ucsc.tables.mod.txt"
	exit
fi

while IFS='' read -r table || [[ -n "$table" ]]; do
	f=$destRegs/$table.bed
	if [ ! -s $f ] ; then
		echo "fetching coordinates from $table"
		`mysql --user=genome --host=genome-mysql.cse.ucsc.edu $genome -A -sre "SELECT chrom, chromStart, chromEnd, name from $table" > $f`
		if [ ! -s $f ] ; then
			`mysql --user=genome --host=genome-mysql.cse.ucsc.edu $genome -A -sre "SELECT chrom, txStart, txEnd, name from $table" > $f`
			if [ ! -s $f ] ; then
				`mysql --user=genome --host=genome-mysql.cse.ucsc.edu $genome -A -sre "SELECT chrom, start, end, name from $table" > $f`
				if [ ! -s $f ] ; then
					`mysql --user=genome --host=genome-mysql.cse.ucsc.edu $genome -A -sre "SELECT genoName, genoStart, genoEnd, repName, repClass from $table" > $f`
				fi
			fi
		fi
	fi
done < $destTmp/ucsc.tables.mod.txt 

cd $destRegs
echo "filename	description" > ../index.txt
ls *bed | awk '{print $1,"	",$1}' >> ../index.txt
