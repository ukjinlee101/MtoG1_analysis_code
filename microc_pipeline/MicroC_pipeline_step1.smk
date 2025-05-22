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
rPE = re.compile(".*_R1.fastq.gz")
fileList = list(filter(rPE.match, os.listdir(config["fastqDir"])))
sampleList = [str.replace("_R1.fastq.gz", "") for str in fileList]

print("### Following list of samples will be processed...")
print(sampleList)

###############################################################################
# RULES
rule all:
	input:
		expand("result/log/success/fastqc_{sample}.ok", sample = sampleList),
		expand("result/log/success/alignment_{sample}.ok", sample = sampleList),
		expand("result/log/success/SAMtoPAIRSAM_{sample}.ok", sample = sampleList),
###############################################################################

rule fastqc:
	conda: "MicroC_QC"
	threads: 4
	input:
		R1 = config["fastqDir"] + "/{sample}_R1.fastq.gz",
		R2 = config["fastqDir"] + "/{sample}_R2.fastq.gz",
	output:
		success = "result/log/success/fastqc_{sample}.ok"
	log:
		R1 = "result/log/fastqc/log_{sample}_R1.txt",
		R2 = "result/log/fastqc/log_{sample}_R2.txt"
	params:
		name = "fastqc_{sample}"
	shell:
		"""
		mkdir -p result/fastqc

		fastqc {input.R1} \
			-o result/fastqc \
			-t {threads} \
			&> {log.R1}

		fastqc {input.R2} \
			-o result/fastqc \
			-t {threads} \
			&> {log.R2}

		# Made it?
		touch {output.success}
		"""

rule alignment:
	conda: "MicroC_QC"
	threads: 16
	input: 
		R1 = config["fastqDir"] + "/{sample}_R1.fastq.gz",
		R2 = config["fastqDir"] + "/{sample}_R2.fastq.gz",
	output:
		SAM = "result/SAM/{sample}.sam",
		success = "result/log/success/alignment_{sample}.ok"
	params:
		BWA_DIR = config["BWA_DIR"],
		name = "alignment_{sample}",
	shell:
		"""

		bwa mem -5SP -T0 -t{threads} {params.BWA_DIR}/genome.fa {input.R1} {input.R2} > {output.SAM}

		# Made it?
		touch {output.success}
		"""

rule SAMtoPAIRSAM:
	conda: "MicroC_QC"
	threads: 16
	input:
		SAM = "result/SAM/{sample}.sam",
	output:
		success = "result/log/success/SAMtoPAIRSAM_{sample}.ok",
		PAIRSAM = "result/pairSAM/{sample}.pairsam",
	params:
		name = "SAMtoPAIRSAM_{sample}",
		REF_DIR = config["REF_DIR"],
		TEMP_DIR = "result/temp"
	shell:
		"""
		mkdir -p result/temp

		pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 \
			--nproc-in {threads} --nproc-out {threads} \
			--chroms-path {params.REF_DIR}/genome.fa {input.SAM} | \
		pairtools sort --tmpdir={params.TEMP_DIR} > {output.PAIRSAM}

		# Made it?
		touch {output.success}
		"""