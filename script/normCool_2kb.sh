mamba activate hicexplorer

sample="G1DMSO_pooled"
sample="G1dTAG_pooled"
sample="G1A485_pooled"
sample="EpiG1DMSO_pooled"
sample="EpiG1dTAG_pooled"
sample="GSE178982_AsyncUT_pooled"
sample="GSE178982_AsyncAID_pooled"

workDir="/athena/apostoloulab/scratch/ukl4001/data/cool_norm_pooled"
cooler coarsen -k 2 -o ${workDir}/${sample}_2000bp_KR.cool ${workDir}/${sample}_1000bp_KR.cool
cooler balance --force -p 32 ${workDir}/${sample}_2000bp_KR.cool
cooltools expected-cis -p 32 -o ${workDir}/${sample}_2000bp_KR_exp.tsv ${workDir}/${sample}_2000bp_KR.cool