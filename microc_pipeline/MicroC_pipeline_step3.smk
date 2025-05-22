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
sampleList = ["G1DMSO",
			  "G1dTAG",
			  "G1A485"]
sampleBatchList = ["G1DMSO_B1", "G1DMSO_B2",
				   "G1dTAG_B1", "G1dTAG_B3",
				   "G1A485_B2", "G1A485_B3"]
sampleRepList = ["G1DMSO_B1R1", "G1DMSO_B1R2", "G1DMSO_B1R3",
				"G1DMSO_B2R1", "G1DMSO_B2R2", "G1DMSO_B2R3",
				"G1A485_B2R1", "G1A485_B2R2", "G1A485_B2R3",
				"G1A485_B3R1", "G1A485_B3R2", "G1A485_B3R3",
				"G1dTAG_B1R1", "G1dTAG_B1R2", "G1dTAG_B1R3",
				"G1dTAG_B3R1", "G1dTAG_B3R2", "G1dTAG_B3R3"]
resList = [1000, 2500, 5000, 10000, 25000, 100000]

###############################################################################
# RULES
rule all:
	input:
		expand("result/log/success/makeHic_{sampleRep}.ok", sampleRep = sampleRepList),
		expand("result/log/success/hictomcool_{sampleRep}.ok", sampleRep = sampleRepList),

		expand("result/log/success/mergePairsBatch_{sampleBatch}.ok", sampleBatch = sampleBatchList),
		expand("result/log/success/makeHicBatch_{sampleBatch}.ok", sampleBatch = sampleBatchList),
		expand("result/log/success/hictomcoolBatch_{sampleBatch}.ok", sampleBatch= sampleBatchList),

		expand("result/log/success/mergePairs_{sample}.ok", sample = sampleList),
		expand("result/log/success/makeHicPooled_{sample}.ok", sample = sampleList),
		expand("result/log/success/hictomcoolPooled_{sample}.ok", sample = sampleList),
###############################################################################

rule makeHic:
	conda: "MicroC_QC"
	threads: 16
	input:
		PAIRS = "result/mappedPair/{sampleRep}.pairs",
	output:
		success = "result/log/success/makeHic_{sampleRep}.ok",
		hic = "result/hic/{sampleRep}_allRes.hic",
	shell:
		"""
		java -Xmx256g -Djava.awt.headless=true -jar /athena/apostoloulab/scratch/ukl4001/001_software/juicer_tools_1.22.01.jar pre \
		--threads {threads} -r 1000,2000,2500,5000,10000,25000,50000,100000,250000,500000,1000000,3000000 \
		{input.PAIRS} {output.hic} /athena/apostoloulab/scratch/collab/Genome/Mus_Musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.sizes

		# Made it?
		touch {output.success}
		"""

rule hictomcool:
	conda: "hicexplorer"
	threads: 16
	input:
		hic = "result/hic/{sampleRep}_allRes.hic"
	output:
		success = "result/log/success/hictomcool_{sampleRep}.ok",
 		mcool = "result/mcool/{sampleRep}_allRes.mcool"
	shell:
		"""
		hicConvertFormat -m {input.hic} --inputFormat hic -o {output.mcool} --outputFormat cool

		# Made it?
		touch {output.success}
		"""

rule mergePairsBatch:
	conda: "MicroC_QC"
	threads: 16
	output:
		pairs = "result/mappedPair_batch/{sampleBatch}.pairs",
		success = "result/log/success/mergePairsBatch_{sampleBatch}.ok",
		log = "result/log/mergePairs/log_{sampleBatch}.txt"
	params:
		name = "{sampleBatch}"
	shell:
		"""
		ls result/mappedPair/{params.name}*.pairs > {output.log}
		
		pairtools merge \
			--nproc {threads} \
			--output {output.pairs} \
			result/mappedPair/{params.name}*.pairs

		# Made it?
		touch {output.success}
		"""

rule makeHicBatch:
	conda: "MicroC_QC"
	threads: 16
	input:
		PAIRS = "result/mappedPair_batch/{sampleBatch}.pairs",
	output:
		success = "result/log/success/makeHicBatch_{sampleBatch}.ok",
		hic = "result/hic_batch/{sampleBatch}_allRes.hic",
	shell:
		"""
		java -Xmx256g -Djava.awt.headless=true -jar /athena/apostoloulab/scratch/ukl4001/001_software/juicer_tools_1.22.01.jar pre \
		--threads {threads} -r 1000,2000,2500,5000,10000,25000,50000,100000,250000,500000,1000000,3000000 \
		{input.PAIRS} {output.hic} /athena/apostoloulab/scratch/collab/Genome/Mus_Musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.sizes

		# Made it?
		touch {output.success}
		"""

rule hictomcoolBatch:
	conda: "hicexplorer"
	threads: 16
	input:
		hic = "result/hic_batch/{sampleBatch}_allRes.hic"
	output:
		success = "result/log/success/hictomcoolBatch_{sampleBatch}.ok",
 		mcool = "result/mcool_batch/{sampleBatch}_allRes.mcool"
	shell:
		"""
		hicConvertFormat -m {input.hic} --inputFormat hic -o {output.mcool} --outputFormat cool

		# Made it?
		touch {output.success}
		"""

rule mergePairs:
	conda: "MicroC_QC"
	threads: 16
	output:
		pairs = "result/mappedPair_pooled/{sample}_pooled.pairs",
		success = "result/log/success/mergePairs_{sample}.ok",
		log = "result/log/mergePairs/log_{sample}.txt"
	params:
		name = "{sample}"
	shell:
		"""
		ls result/mappedPair/{params.name}*.pairs > {output.log}
		
		pairtools merge \
			--nproc {threads} \
			--output {output.pairs} \
			result/mappedPair/{params.name}*.pairs


		# Made it?
		touch {output.success}
		"""

rule makeHicPooled:
	conda: "MicroC_QC"
	threads: 16
	input:
		PAIRS = "result/mappedPair_pooled/{sample}_pooled.pairs",
	output:
		success = "result/log/success/makeHicPooled_{sample}.ok",
		hic = "result/hic_pooled/{sample}_pooled_allRes.hic",
	shell:
		"""
		java -Xmx512g -Djava.awt.headless=true -jar /athena/apostoloulab/scratch/ukl4001/001_software/juicer_tools_1.22.01.jar pre \
		--threads {threads} -r 1000,2000,2500,5000,10000,25000,50000,100000,250000,500000,1000000,3000000 \
		{input.PAIRS} {output.hic} /athena/apostoloulab/scratch/collab/Genome/Mus_Musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.sizes

		# Made it?
		touch {output.success}
		"""

rule hictomcoolPooled:
	conda: "hicexplorer"
	threads: 16
	input:
		hic = "result/hic_pooled/{sample}_pooled_allRes.hic"
	output:
		success = "result/log/success/hictomcoolPooled_{sample}.ok",
 		mcool = "result/mcool_pooled/{sample}_pooled_allRes.mcool"
	shell:
		"""
		hicConvertFormat -m {input.hic} --inputFormat hic -o {output.mcool} --outputFormat cool

		# Made it?
		touch {output.success}
		"""