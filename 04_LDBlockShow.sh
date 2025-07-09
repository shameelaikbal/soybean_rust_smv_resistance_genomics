#LDblockShow
LDBlockShow-1.40/bin/LDBlockShow -InVCF region_around_GWASNP.vcf -SeleVar 3 -BlockType 2 -Region glyma.Wm82.gnm2.Gm18:56098721-56298006 \
-InGWAS region_farmcpu_cleaned.txt -OutPng -OutPut region_new_200kb_chr18_56100116_ingwas

 LDBlockShow-1.40 % ./bin/ShowLDSVG \     
    -InPreFix region_new_200kb_chr18_56100116_ingwas \
    -InGWAS region_farmcpu_cleaned.txt \
    -OutPut LD_white_to_purple \
    -crBegin 255,255,255 \
    -crMiddle 255,127,0 \
    -crEnd 128,0,128 \
    -TopSite glyma.Wm82.gnm2.Gm18:56100116 \
    -SpeSNPName Spe.snp \    
    -ShowGWASSpeSNP  \
    -OutPng
