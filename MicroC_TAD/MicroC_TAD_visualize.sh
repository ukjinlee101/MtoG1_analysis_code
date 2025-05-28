mamba activate jupyterlab_cooltools

COOL='/Volumes/UKJIN_SSD/data/cool_norm_pooled/G1DMSO_pooled_25000bp_KR.cool'
TAD='/Volumes/UKJIN_SSD/data/TAD/TAD_25kb_125kb_otsu_G1DMSO_pooled.bedpe'

coolpup.py --features_format bedpe \
            --ignore_diags 2 \
            --local \
            --rescale \
            --rescale_size 99 \
            --rescale_flank 1 \
            --clr_weight_name weight \
            --nproc 10 \
            --seed 123 \
            "${COOL}" "${TAD}"
            


pileup_df = coolpup.pileup(clr = clr,
                        features = tads,
                        features_format = 'bedpe',
                        clr_weight_name = 'weight',
                        min_diag = 2,
                        local = True,
                        rescale = True,
                        rescale_size = 99,
                        rescale_flank = 1,
                        nproc = 10, seed = 123)