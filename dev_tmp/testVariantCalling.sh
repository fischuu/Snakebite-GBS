/scripts/GBS-SNP-CROP-7.pl -in {params.input} \
                           -out {params.out} \
                            -p {params.p} \
                            -mnHoDepth0 {params.mnHoDepth0} \
                            -mnHoDepth1 {params.mnHoDepth1} \
                            -mnHetDepth {params.mnHetDepth} \
                            -altStrength {params.altStrength} \
                            -mnAlleleRatio {params.mnAlleleRatio} \
                            -mnCall {params.mnCall} \
                            -mnAvgDepth {params.mnAvgDepth} \
                            -mxAvgDepth {params.mxAvgDepth}

perl /users/fischerd/git/Snakebite-GBS/scripts/GBS-SNP-CROP-8_1.pl -in GSC.GenoMatrix.txt \
                                                            -out test.vcf \
                                                            -b barcodesID.txt \
                                                            -formats vcf