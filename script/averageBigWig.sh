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

