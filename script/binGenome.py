#!/usr/bin/env python3

import gzip
import os
import sys

if len(sys.argv) != 4:
    sys.exit(f"Usage: {sys.argv[0]} <genomeSize> <binSize> <outFile>")

genomeSizeFile, binSize, outFile = sys.argv[1], int(sys.argv[2]), sys.argv[3]

# FUNCTION
def generate_bins(genomeSizeFile, binSize, outputFile):
	bins = []
	with open(genomeSizeFile, 'r') as f:
		for line in f:
			chrom, size = line.strip().split()
			size = int(size)
			for start in range (0, size, binSize):
				end = min(start + binSize, size)
				binName = start
				bins.append((chrom, start, end, binName))
	with gzip.open(outputFile, 'wt') as out:
		for bin in bins:
			out.write(f"{bin[0]}\t{bin[1]}\t{bin[2]}\t{bin[3]}\n")

# MAIN
generate_bins(genomeSizeFile, binSize, outFile)