# INSTALLING GENOMEDISCO
# mamba create -n genomedisco_python2
# mamba activate genomedisco_python2
# conda install python=2.7 
# cd ~/Local_Bioinformatics/001_software
# git clone http://github.com/kundajelab/genomedisco
# mamba install pip --force-reinstall
# pip install --editable genomedisco
# python -m pip install scipy
# python -m pip install -U scikit-learn
# python -m pip install psutil

# MAKING BED FILE OF ALL BINS USED IN THE ANALYSIS
SCRDIR="../script"
REFDIR="../reference/mm10"
"$SCRDIR"/binGenome.py  ${REFDIR}/mm10.chrom.sizes.canonical 50000 ${REFDIR}/mm10.bin.50kb.bed.gz

# MAKE BINNED MATRIX FROM HIC
./binHiC.sh

# MAKE METADATA
./makeMetaData.sh

# RUN
sbatch runSnakemake.sh MicroC_GenomeDISCO_pipeline.smk