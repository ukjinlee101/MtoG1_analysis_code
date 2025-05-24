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
		#expand("result/log/success/normCool_{sampleRep}_{res}.ok", sampleRep = sampleRepList, res = resList),
		expand("result/log/success/normCoolBatch_{sampleBatch}_{res}.ok", sampleBatch = sampleBatchList, res = resList),
		expand("result/log/success/normCoolPooled_{sample}_{res}.ok", sample = sampleList, res = resList),
###############################################################################

rule normCool:
	conda: "hicexplorer"
	threads: 16
	input:
		mcool = "result/mcool/{sampleRep}_allRes.mcool"
	output:
		success = "result/log/success/normCool_{sampleRep}_{res}.ok",
		mcool = "result/mcool_norm/{sampleRep}_{res}bp_KR.mcool",
		exp = "result/mcool_norm/{sampleRep}_{res}bp_KR_exp.tsv"
	params:
		res = "{res}"
	shell:
		"""

		hicConvertFormat -m {input.mcool}:://resolutions/{params.res} --inputFormat cool \
			--outputFormat cool -o {output.mcool} --correction_name KR
		cooltools expected-cis -p 16 -o {output.exp} {output.mcool}

		# Made it?
		touch {output.success}
		"""

rule normCoolBatch:
	conda: "hicexplorer"
	threads: 16
	input:
		mcool = "result/mcool_batch/{sampleBatch}_allRes.mcool"
	output:
		success = "result/log/success/normCoolBatch_{sampleBatch}_{res}.ok",
		mcool = "result/mcool_norm_batch/{sampleBatch}_{res}bp_KR.mcool",
		exp = "result/mcool_norm_batch/{sampleBatch}_{res}bp_KR_exp.tsv"
	params:
		res = "{res}"
	shell:
		"""

		hicConvertFormat -m {input.mcool}:://resolutions/{params.res} --inputFormat cool \
			--outputFormat cool -o {output.mcool} --correction_name KR
		cooltools expected-cis -p 16 -o {output.exp} {output.mcool}

		# Made it?
		touch {output.success}
		"""

rule normCoolPooled:
	conda: "hicexplorer"
	threads: 16
	input:
		mcool = "result/mcool_pooled/{sample}_pooled_allRes.mcool"
	output:
		success = "result/log/success/normCoolPooled_{sample}_{res}.ok",
		mcool = "result/mcool_norm_pooled/{sample}_pooled_{res}bp_KR.mcool",
		exp = "result/mcool_norm_pooled/{sample}_pooled_{res}bp_KR_exp.tsv"
	params:
		res = "{res}"
	shell:
		"""

		hicConvertFormat -m {input.mcool}:://resolutions/{params.res} --inputFormat cool \
			--outputFormat cool -o {output.mcool} --correction_name KR
		cooltools expected-cis -p 16 -o {output.exp} {output.mcool}

		# Made it?
		touch {output.success}
		"""