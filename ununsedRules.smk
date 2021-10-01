if config["mockreference"] == "":
    pass
else:
    rule ParseMpileup_createCountFiles_existingMockreference:
        """
        Parse mpileup outputs and create count/ref files for existing mock reference genome
        """
        input:
            "%s/MPILEUP/mpileup_existingMock/{samples}.mpileup" % (config["project-folder"])
        output:
            co="%s/MPILEUP/mpileup_existingMock/{samples}.count.txt" % (config["project-folder"]),
            ref="%s/MPILEUP/mpileup_existingMock/{samples}.ref.txt" % (config["project-folder"])
        params:
            p=config["params"]["step6"]["p"],
            wd="%s/MPILEUP/mpileup_existingMock" % config["project-folder"],
            pipefolder=config["pipeline-folder"]
        log:
            "%s/logs/Perl/ParseMpileup_createCountFiles_exmockreference.{samples}.log" % (config["project-folder"])
        benchmark:
            "%s/benchmark/Perl/ParseMpileup_createCountFiles_exmockreference.{samples}.benchmark.tsv" % (config["project-folder"])
        threads: 1
        singularity: config["singularity"]["gbs"]
        shell:"""
            cd {params.wd}
            perl {params.pipefolder}/scripts/GBS-SNP-CROP-6_1.pl -b {input} -p {params.p} &> {log}
      	"""

if config["mockreference"] == "":
    pass
else:
    rule create_verticalRef_existingMockreference:
        """
        Merge the separate ref files for existing mock reference genome
        """
        input:
            co=expand("%s/MPILEUP/mpileup_existingMock/{samples}.count.txt" % (config["project-folder"]), samples=samples),
            ref=expand("%s/MPILEUP/mpileup_existingMock/{samples}.ref.txt" % (config["project-folder"]), samples=samples),
        output:
           vref="%s/MPILEUP/mpileup_existingMock/VerticalRefPos.txt" % (config["project-folder"]),
           co="%s/MPILEUP/mpileup_existingMock/CountFileList.txt" % (config["project-folder"])
        log:
           "%s/logs/BASH/create_verticalRef_exmockreference.log" % (config["project-folder"])
        benchmark:
            "%s/benchmark/BASH/create_verticalRef_exmockreference.benchmark.tsv" % (config["project-folder"])
        params:
            wd="%s/MPILEUP/mpileup_existingMock" % config["project-folder"],
        threads: 1
        shell:"""
            sort -u -m -k2n {input.ref} -o {params.wd}/VerticalRefPos.tmp;
            sort -k2n {params.wd}/VerticalRefPos.tmp | uniq > {output.vref};
           # ls {params.wd}/*count.txt | xargs -n1 basename > {output.co}
            ls {params.wd}/*count.txt > {output.co}
        """
        
if config["mockreference"] == "":
    pass
else:
    checkpoint cut_verticalRef_existingMockreference:
        """
        Divide the input verticalRef-file for parallel processing
        """
        input:
           "%s/MPILEUP/mpileup_existingMock/VerticalRefPos.txt" % (config["project-folder"])
        output:
            directory("%s/MPILEUP/mpileup_existingMock/VerticalRefPos/" % (config["project-folder"]))
        params:
            out="%s/MPILEUP/mpileup_existingMock/VerticalRefPos/VerticalRefPos." % (config["project-folder"]),
            split=1000000
        shell:"""
            mkdir -p {output}
            split -l {params.split} --numeric-suffixes {input} {params.out}
        """
        
if config["mockreference"] == "":
    pass
else:
    rule create_MasterMatrix_parallel_existingMockreference:
        """
        Process verticalRef parallel to create MasterMatrix
        """
        input:
            verRef="%s/MPILEUP/mpileup_existingMock/VerticalRefPos/VerticalRefPos.{i}" % (config["project-folder"]),
            counts="%s/MPILEUP/mpileup_existingMock/CountFileList.txt" % (config["project-folder"])
        output:
            "%s/MPILEUP/mpileup_existingMock/VerticalRefPos/GSC.MasterMatrix_{i}.tsv" % (config["project-folder"])
        log:
            "%s/logs/BASH/create_MasterMatrix_parallel_exmockreference_{i}.log" % (config["project-folder"])
        benchmark:
            "%s/benchmark/BASH/create_MasterMatrix_parallel_exmockreference.benchmark_{i}.tsv" % (config["project-folder"])
        params:
            p=config["params"]["step6"]["p"],
            wd="%s/MPILEUP/mpileup_existingMock/VerticalRefPos" % config["project-folder"],
            pipefolder=config["pipeline-folder"]
        conda:"envs/gbs.yaml"
        singularity: config["singularity"]["gbs"]
        shell:"""
            cd {params.wd}
            perl {params.pipefolder}/scripts/GBS-SNP-CROP-6_2.pl -in {input.verRef} -count {input.counts} -out {output} -p {params.p} &> {log}
         """
         
if config["mockreference"] == "":
    pass
else:
    def aggregate_inputMasterMatrix_existingMockreference(wildcards):
        """
        Aggregate the input object for the final MasterMatrix
        """
        checkpoint_outputMM = checkpoints.cut_verticalRef_existingMockreference.get(**wildcards).output[0]
        return expand("%s/MPILEUP/mpileup_existingMock/VerticalRefPos/GSC.MasterMatrix_{i}.tsv" % (config["project-folder"]),
                      i=glob_wildcards(os.path.join(checkpoint_outputMM, "VerticalRefPos.{i}")).i)      
                      

if config["mockreference"] == "":
    pass
else:
    rule aggregate_MasterMatrix_existingMockreference:
        input:
            aggregate_inputMasterMatrix_existingMockreference
        output:
            "%s/MPILEUP/mpileup_existingMock/GSC.MasterMatrix.txt" % (config["project-folder"])
        shell:"""
            cat {input} > {output}
        """   
        
if config["mockreference"] == "":
    pass
else:
    rule FilterVariants_existingMockreference:
        """
        Filter variants and call genotypes using the existing mock reference genome.
        """
        input:
            "%s/MPILEUP/mpileup_existingMock/GSC.MasterMatrix.txt" % (config["project-folder"])
        output:
            filtered="%s/MPILEUP/mpileup_existingMock/variants/GSC.GenoMatrix.txt" % (config["project-folder"]),
            unfiltered="%s/MPILEUP/mpileup_existingMock/variants/GSC.GenoMatrix.unfiltered.txt" % (config["project-folder"])
        params:
            input=config["params"]["step7"]["input"],
            out=config["params"]["step7"]["out"],
            unfiltered=config["params"]["step7"]["unfiltered"],
            p=config["params"]["step7"]["p"],
            mnHoDepth0=config["params"]["step7"]["mnHoDepth0"],
            mnHoDepth1=config["params"]["step7"]["mnHoDepth1"],
            mnHetDepth=config["params"]["step7"]["mnHetDepth"],
            altStrength=config["params"]["step7"]["altStrength"],
            mnAlleleRatio=config["params"]["step7"]["mnAlleleRatio"],
            mnCall=config["params"]["step7"]["mnCall"],
            mnAvgDepth=config["params"]["step7"]["mnAvgDepth"],
            mxAvgDepth=config["params"]["step7"]["mxAvgDepth"],
            wd="%s/MPILEUP/mpileup_existingMock" % config["project-folder"],
            pipefolder=config["pipeline-folder"]
        log:
            "%s/logs/Perl/FilterVariants_exMockreference.log" % (config["project-folder"])
        benchmark:
            "%s/benchmark/Perl/FilterVariants_exMockreference.benchmark.tsv" % (config["project-folder"])
        singularity: config["singularity"]["gbs"]
        shell:"""
            #******PARAMETERS*****
            # -in 	Discovery master matrix input file. The output from last step 	File 	GSC.MasterMatrix.txt
            # -out 	Genotyping matrix for the population 	Output file 	GSC.GenoMatrix.txt
            # -p 	Identify SNPs only (snp) or SNPs + indels (indel) 	String 	--
            # -mnHoDepth0 	Minimum depth required for calling a homozygote when the alternative allele depth = 0 	Numeric 	5
            # -mnHoDepth1 	Minimum depth required for calling a homozygote when the alternative allele depth = 1 	Numeric 	20
            # -mnHetDepth 	Minimum depth required for each allele when calling a heterozygote 	Numeric 	3
            # -altStrength 	Across the population for a given putative bi-allelic variant, this alternate allele strength parameter is the minimum proportion of non-primary allele reads that are the secondary allele 	Numeric 	0.8
            # -mnAlleleRatio 	Minimum required ratio of less frequent allele depth to more frequent allele depth 	Numeric 	0.25
            # -mnCall 	Minimum acceptable proportion of genotyped individuals to retain a variant 	Numeric 	0.75
            # -mnAvgDepth 	Minimum average depth of an acceptable variant 	Numeric 	3
            # -mxAvgDepth 	Maximum average depth of an acceptable variant 	Numeric 	200
    
            cd {params.wd}
            perl {params.pipefolder}/scripts/GBS-SNP-CROP-7.pl -in {params.input} -out {params.out} -p {params.p} -mnHoDepth0 {params.mnHoDepth0} -mnHoDepth1 {params.mnHoDepth1} -mnHetDepth {params.mnHetDepth} -altStrength {params.altStrength} -mnAlleleRatio {params.mnAlleleRatio} -mnCall {params.mnCall} -mnAvgDepth {params.mnAvgDepth} -mxAvgDepth {params.mxAvgDepth} &> {log}
            perl {params.pipefolder}/scripts/GBS-SNP-CROP-7.pl -in {params.input} -out {params.unfiltered} -p {params.p} -mnHoDepth0 0 -mnHoDepth1 0 -mnHetDepth 0 -altStrength 0 -mnAlleleRatio 0 -mnCall 0 -mnAvgDepth 0 -mxAvgDepth 1000000 &> {log}
      	"""
       
if config["mockreference"] == "":
    pass
else:
    rule CreateVCF_existingMockreference:
        """
        Create the VCF output of the pipeline for the existing mock reference genome.
        """
        input:
            filtered="%s/MPILEUP/mpileup_existingMock/variants/GSC.GenoMatrix.txt" % (config["project-folder"]),
            unfiltered="%s/MPILEUP/mpileup_existingMock/variants/GSC.GenoMatrix.unfiltered.txt" % (config["project-folder"])
        output:
            "%s/VCF/FinalSetVariants_existingMock.vcf" % (config["project-folder"])
        params:
            barcodes=config["barcodes"],
            genoMatrix=config["params"]["step8"]["genoMatrix"],
            out=config["params"]["step8"]["out"],
            unfilteredOut=config["params"]["step8"]["unfiltered"],
            format=config["params"]["step8"]["format"],
            wd="%s/MPILEUP/mpileup_existingMock" % config["project-folder"],
            outdir="%s/VCF" % config["project-folder"],
            pipeout="%s/MPILEUP/mpileup_existingMock/GSC.vcf" % (config["project-folder"]),
            pipefolder=config["pipeline-folder"]
        log:
            "%s/logs/Perl/CreateVCF_exMockreference.log" % (config["project-folder"])
        benchmark:
            "%s/benchmark/Perl/CreateVCF_exMockreference.benchmark.tsv" % (config["project-folder"])
        singularity: config["singularity"]["gbs"]
        shell:"""
            #******PARAMETERS*****
            # -in 	Final genotyping matrix input file 	File 	GSC.GenoMatrix.txt
            # -out 	Output label name without extension 	String 	GSC
            # -b 	BarcodeID file name See Appendix A 	File 	--
            # -formats 	The name(s) (R, Tassel, Plink, vcf, H) for which the GBS-SNP-CROP final genotyping matrix should be converted. If more than one format is desired, the names should be separated by commas without any space, as in the examples shown below 	String 	R,T,P,V,H
    
            cd {params.wd}
            perl {params.pipefolder}/scripts/GBS-SNP-CROP-8.pl -in {input.filtered} -out {params.out} -b {params.barcodes} -formats {params.format} &> {log}
            perl {params.pipefolder}/scripts/GBS-SNP-CROP-8.pl -in {input.unfiltered} -out {params.unfilteredOut} -b {params.barcodes} -formats {params.format} &> {log}
            cp {params.pipeout} {output}
            cd ../..
      	"""