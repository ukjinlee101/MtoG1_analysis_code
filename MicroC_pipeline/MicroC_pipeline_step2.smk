###############################################################################
### Micro-C pipeline for PE dataset                                         ###
###############################################################################

###############################################################################
# REFERENCE
# https://micro-c.readthedocs.io/en/latest/before_you_begin.html
###############################################################################

### Importing python packages
import os
import pandas as pd
import re

### Configuration
configfile: "config/config.json"

### Extracting sample names
sampleList = ["G1DMSO_B1R1", "G1DMSO_B1R2", "G1DMSO_B1R3",
			  "G1DMSO_B2R1", "G1DMSO_B2R2", "G1DMSO_B2R3",
			  "G1dTAG_B1R1", "G1dTAG_B1R2", "G1dTAG_B1R3",
			  "G1dTAG_B3R1", "G1dTAG_B3R2", "G1dTAG_B3R3",
			  "G1A485_B2R1", "G1A485_B2R2", "G1A485_B2R3",
			  "G1A485_B3R1", "G1A485_B3R2", "G1A485_B3R3"]
print("### Following list of samples will be processed...")
print(sampleList)

###############################################################################
# RULES
rule all:
	input:
		expand("result/log/success/mergePairsam_{sample}.ok", sample = sampleList),
		expand("result/log/success/dedup_{sample}.ok", sample = sampleList),
		expand("result/log/success/dedupQC_{sample}.ok", sample = sampleList),
		expand("result/log/success/classifyDup_{sample}.ok", sample = sampleList),
###############################################################################

rule mergePairsam:
	conda: "MicroC_QC"
	threads: 16
	output:
		pairsam = "result/pairSAM_pooled/{sample}_pooled.pairsam",
		success = "result/log/success/mergePairsam_{sample}.ok",
		log = "result/log/mergePairsam/log_{sample}.txt"
	params:
		name = "{sample}"
	shell:
		"""
		# Log input files
		ls result/pairSAM/{params.name}*.pairsam > {output.log}

		pairtools merge \
			--nproc {threads} \
			--output {output.pairsam} \
			result/pairSAM/{params.name}*.pairsam

		# Made it?
		touch {output.success}
		"""

rule dedup:
	conda: "MicroC_QC"
	threads: 16
	input:
		PAIRSAM = "result/pairSAM_pooled/{sample}_pooled.pairsam",
	output:
		success = "result/log/success/dedup_{sample}.ok",
		STATS = "result/stats/stats_{sample}.txt",
		DEDUPPAIRSAM = temp("result/temp/{sample}_dedup.pairsam"),
		MAPPEDPAIR = "result/mappedPair/{sample}.pairs",
		MAPPEDUNSBAM = temp("result/temp/{sample}_mappeduns.bam"),
		MAPPEDBAM = "result/mappedBAM/{sample}.bam",
		MAPPEDBAMIndex = "result/mappedBAM/{sample}.bam.bai",
	params:
		name = "dedup_{sample}",
	shell:
		"""
		
		mkdir -p result/temp

		pairtools dedup --nproc-in {threads} --nproc-out {threads} --mark-dups \
			--output-stats {output.STATS} --output {output.DEDUPPAIRSAM} {input.PAIRSAM}


		pairtools split --nproc-in {threads} --nproc-out {threads} --output-pairs {output.MAPPEDPAIR} --output-sam {output.MAPPEDUNSBAM} {output.DEDUPPAIRSAM}

		samtools sort -@{threads} -o {output.MAPPEDBAM} {output.MAPPEDUNSBAM}
		samtools index {output.MAPPEDBAM}

		# Made it?
		touch {output.success}
		"""

rule dedupQC:
	conda: "MicroC_QC"
	threads: 4
	input:
		STATS = "result/stats/stats_{sample}.txt",
	output:
		success = "result/log/success/dedupQC_{sample}.ok",
		STATS_SUM = "result/stats/summary_{sample}.txt",
	shell:
		"""
		python3 get_qc.py -p {input.STATS} > {output.STATS_SUM}

		# Made it?
		touch {output.success}
		"""

rule classifyDup:
	conda: "MicroC_QC"
	threads: 16
	input:
		PAIRSAM = "result/pairSAM_pooled/{sample}_pooled.pairsam",
	output:
		success = "result/log/success/classifyDup_{sample}.ok",
		NONDEDUPPAIR = temp("result/temp/{sample}_nondedup.pairs"),
		NONDEDUPBAM = temp("result/temp/{sample}_nondedup.bam"),
		namecollate = temp("result/temp/{sample}_namecollate.bam"),
		fixmate = temp("result/temp/{sample}_fixmate.bam"),
		positionsort = temp("result/temp/{sample}_positionsort.bam"),
		NONDEDUPSTATSUM = "result/stats/summary_{sample}_duplicates.txt",
		MARKDUPBAM = temp("result/temp/{sample}_markdup.bam"),
	params:
		name = "classifyDup_{sample}",
	shell:
		"""
		
		pairtools split --nproc-in {threads} --nproc-out {threads} \
			--output-pairs {output.NONDEDUPPAIR} \
			--output-sam {output.NONDEDUPBAM} \
			{input.PAIRSAM}

		samtools collate -@ {threads} -o {output.namecollate} {output.NONDEDUPBAM}
		samtools fixmate -@ {threads} -m {output.namecollate} {output.fixmate}
		samtools sort -@ {threads} -o {output.positionsort} {output.fixmate}
		samtools markdup -@ {threads} -f {output.NONDEDUPSTATSUM} -S -d 2500 -m s --include-fails {output.positionsort} {output.MARKDUPBAM}

		# Made it?
		touch {output.success}
		"""

		