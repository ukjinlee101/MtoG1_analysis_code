mamba activate bioinfo

wiggletools mean H3K4me3_ESC_Rep1_depthNorm.bw H3K4me3_ESC_Rep2_depthNorm.bw H3K4me3_ESC_Rep3_depthNorm.bw > temp.wig
wigToBigWig temp.wig ../mm10/mm10.chrom.sizes H3K4me3_ESC_avg_depthNorm.bw
rm temp.wig

wiggletools mean H3K27ac_ESC_Rep1_depthNorm.bw H3K27ac_ESC_Rep2_depthNorm.bw H3K27ac_ESC_Rep3_GSM2438476_depthNorm.bw > temp.wig
wigToBigWig temp.wig ../mm10/mm10.chrom.sizes H3K27ac_ESC_avg_depthNorm.bw
rm temp.wig

wiggletools mean H3K4me3_EpiLC_Rep1_depthNorm.bw H3K4me3_EpiLC_Rep2_depthNorm.bw > temp.wig
wigToBigWig temp.wig ../../mm10/mm10.chrom.sizes H3K4me3_EpiLC_avg_depthNorm.bw
rm temp.wig

wiggletools mean H3K27ac_EpiLC_Rep1_depthNorm.bw H3K27ac_EpiLC_Rep2_depthNorm.bw > temp.wig
wigToBigWig temp.wig ../../mm10/mm10.chrom.sizes H3K27ac_EpiLC_avg_depthNorm.bw
rm temp.wig


wiggletools mean RAD21_ESC-DMSO_Rep1_spikeNorm_inputNorm.bw RAD21_ESC-DMSO_Rep2_spikeNorm_inputNorm.bw > temp.wig
wigToBigWig temp.wig ../../mm10/mm10.chrom.sizes RAD21_ESC-DMSO_avg_spikeNorm_inputNorm.bw
rm temp.wig
wiggletools mean RAD21_EpiLC-DMSO_Rep1_spikeNorm_inputNorm.bw RAD21_EpiLC-DMSO_Rep2_spikeNorm_inputNorm.bw > temp.wig
wigToBigWig temp.wig ../../mm10/mm10.chrom.sizes RAD21_EpiLC-DMSO_avg_spikeNorm_inputNorm.bw
rm temp.wig
