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
sampleList = ["GSE178982_AsyncUT_pooled", "GSE178982_AsyncAID_pooled", "G1DMSO_pooled", "G1dTAG_pooled", "G1A485_pooled", "EpiG1DMSO_pooled", "EpiG1dTAG_pooled"]
resList = [5000, 2000, 1000]
###############################################################################
# RULES
rule all:
	input:
		expand("/athena/apostoloulab/scratch/ukl4001/data/loop_Hansen/{sample}_union_{res}bp_pu100pz100.tsv", sample = sampleList, res = resList),


###############################################################################
rule calChromoScore:
	conda: "chromosight"
	threads: 16
	input:
		normCool = "/athena/apostoloulab/scratch/ukl4001/data/cool_norm_pooled/{sample}_{res}bp_KR.cool",
		loop = "/athena/apostoloulab/scratch/ukl4001/data/loop_Hansen/Hansen_union_allRes_postprocessed_{res}bp.bedpe",
	output:
		score = "/athena/apostoloulab/scratch/ukl4001/data/loop_Hansen/{sample}_union_{res}bp_pu100pz100.tsv"
	params:
		name = "/athena/apostoloulab/scratch/ukl4001/data/loop_Hansen/{sample}_union_{res}bp_pu100pz100"
	shell:
		"""
        chromosight quantify {input.loop} \
                        {input.normCool} \
                        {params.name} \
                        --pattern=loops \
                        --threads=16 \
                        --perc-undetected=100 \
                        --perc-zero=100
		"""
