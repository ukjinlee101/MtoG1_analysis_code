#! /bin/bash -l
 
#SBATCH --partition=scu-cpu   # cluster-specific
#SBATCH --time=24:00:00   # HH/MM/SS
#SBATCH --ntasks=1            # Number of tasks (processes)
#SBATCH --cpus-per-task=32 
#SBATCH --mem=256Gb   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=hpc_alert_ukl4001@outlook.com
#SBATCH --mail-type=ALL
#SBATCH -e snakeSlurm.err.%J
#SBATCH -o snakeSlurm.out.%J

source ~/.bashrc

PREFIX="G1dTAG_pooled"
REFDIR='/athena/apostoloulab/scratch/ukl4001/MtoG1_analysis_code/reference/mm10'
RESULTDIR='/athena/apostoloulab/scratch/ukl4001/data/cscore'
WINDOW=${REFDIR}/mm10.bin.10kb.bed
INPUTDIR='/athena/apostoloulab/scratch/ukl4001/temp/hicSummary'
SESSION=32
MINDIST=1000000

mkdir -p ${RESULTDIR}

CscoreTool1.1 ${WINDOW} ${INPUTDIR}/${PREFIX}.summary ${RESULTDIR}/${PREFIX}_100kb ${SESSION} ${MINDIST}
