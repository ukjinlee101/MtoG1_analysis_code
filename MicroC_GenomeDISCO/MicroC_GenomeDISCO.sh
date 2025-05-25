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
"$SCRDIR"/binGenome.py  ${REFDIR}/mm10.chrom.sizes.canonical 100000 ${REFDIR}/mm10.bin.100kb.bed.gz

# MAKE BINNED MATRIX FROM HIC
240623_dumpHic

# Run correlation locally
cd /Volumes/UKJIN_SSD/Genomics_03_Analysis_Working/data_vault_2024summer
genomedisco run_all --metadata_samples source/genomedisco/metadata.samples --metadata_pairs source/genomedisco/metadata/metadata.pairs --bins reference/mm10.bin.10kb.bed.gz --outdir result/genomedisco
genomedisco summary --metadata_samples source/genomedisco/metadata.samples --metadata_pairs source/genomedisco/metadata/metadata.pairs --bins reference/mm10.bin.10kb.bed.gz --outdir result/genomedisco
genomedisco cleanup --outdir result/genomedisco

# Run correlation on cayuga
000_rcloneToBox.sh Box:/Genomics_00_Archive_UkJin/MicroC_240613_16465-merged_G1A485/result/hic ./hic
000_rcloneToBox.sh Box:/Genomics_00_Archive_UkJin/MicroC_231120_14995-15390_RAD21/result/hic ./hic_rad21
sbatch genomedisco_concordance_1.sh
