### Importing python packages
import os
import pandas as pd
import re

### Configuration
configfile: "config/config.json"

### Extracting sample names
#rPE = re.compile(".*\.bed$")
fileList = [filename for filename in os.listdir(config["VPDIR"]) if filename.endswith('.bed') and not filename.startswith('._')]

sampleList = [file.removeprefix("vp.").removesuffix(".bed") for file in fileList]
# sampleList = ["ENSMUSG00000028078.Dclk2"]
expList = ["G1DMSO_pooled", "G1dTAG_pooled", "G1A485_pooled", "EpiG1DMSO_pooled", "EpiG1dTAG_pooled", "GSE178982_AsyncUT_pooled", "GSE178982_AsyncAID_pooled"]
resList = ["5", "10"]

print("### Following list of samples will be processed...")
print(sampleList)
print(expList)
print(resList)

###############################################################################
# RULES
rule all:
	input:
		expand(config["TEMPDIR"] + "/step1.{gene}__{res}kb", gene = sampleList, res = resList),
		expand(config["AIDIR"] + "/allinfo.{gene}__{exp}__{res}kb", exp = expList, gene = sampleList, res = resList),

###############################################################################
rule step:
	conda: "bioinfo"
	threads: 4
	input:
		binned = config["REFDIR"] + "/mm10/mm10.bin.{res}kb.bed",
		vp = config["VPDIR"] + "/vp.{gene}.bed",
		genome = config["GENOMEDIR"] + "/genome.fa.fai",
		mboi = config["REFDIR"] + "/mm10/MboI_mm10.bed"
	output:
		region1 = temp(config["TEMPDIR"] + "/region1.{gene}__{res}kb"),
		step1 = temp(config["TEMPDIR"] + "/step1.{gene}__{res}kb"),
		slop = temp(config["TEMPDIR"] + "/slop.{gene}__{res}kb"),
		step = temp(config["TEMPDIR"] + "/step.{gene}__{res}kb"),
	params:
		res = "{res}"
	shell:
		"""
		# Extract bin with VP (region1)
		intersectBed -a {input.binned} -b {input.vp} -wa -u > {output.region1}
		
		# (ARCHIVED) Expand the region of interest to +- 2Mb (slop)
		# slopBed -i {output.region1} -g {input.genome} -b 2000000 > {output.slop}
		
		# Get the chromosome name from the viewpoint file, find the full size of ghromozome, and create slope for the entire chromosome
		grep -w "$(cut -f1 {output.region1})" {input.genome} | awk '{{print $1 "\\t0\\t" $2}}' > {output.slop}


		# Create (res) kb windows with 0.05*(res) bp steps within regions (step)
		window=$((1000*{params.res}))
		step=$((50*{params.res}))
		bedtools makewindows -b {output.slop} -w ${{window}} -s ${{step}} > {output.step}
		
		# Count number of restriction sites within each window (step1)
		intersectBed -a {output.step} -b {input.mboi} -c > {output.step1}
		"""

rule allInfo:
	conda: "bioinfo"
	threads: 4
	input:
		vp = config["VPDIR"] + "/vp.{gene}.bed",
		region1 = config["TEMPDIR"] + "/region1.{gene}__{res}kb",
		step1 = config["TEMPDIR"] + "/step1.{gene}__{res}kb",
	output:
		pairs = temp(config["TEMPDIR"] + "/pairs.{gene}__{exp}__{res}kb"),
		left = temp(config["TEMPDIR"] + "/left.{gene}__{exp}__{res}kb"),
		right = temp(config["TEMPDIR"] + "/right.{gene}__{exp}__{res}kb"),
		singlereads = temp(config["TEMPDIR"] + "/singlereads.{gene}__{exp}__{res}kb"),
		enrichment = temp(config["TEMPDIR"] + "/enrichment.{gene}__{exp}__{res}kb"),
		allinfo = config["AIDIR"] + "/allinfo.{gene}__{exp}__{res}kb",
	params:
		PAIRDIR = config["PAIRDIR"],
		EXP = "{exp}",
		res = "{res}"
	shell:
		"""
		firstline=$(head -n 1 {input.vp})
		chr=$(echo ${{firstline}} | awk '{{print $1}}')

		pairToBed -a {params.PAIRDIR}/{params.EXP}.${{chr}}.bedpe -b {input.region1} -type "either" > {output.pairs}

		cut -f1,2,3 {output.pairs} > {output.left}
		cut -f4,5,6 {output.pairs} > {output.right}
		cat {output.left} {output.right} > {output.singlereads}
		
		#calculate read coverage for all bins around the region of interest
		intersectBed -a {input.step1} -b {output.singlereads} -c > {output.enrichment}
		
		# add location of region of interest - We will delete this from downstream analysis for plots - it is useful to keep it as a separate column.
		intersectBed -a {output.enrichment} -b {input.region1} -c > {output.allinfo}
		
		"""