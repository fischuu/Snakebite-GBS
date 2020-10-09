     # Extract consecutive bases and put them into a bed file with +-6 offset from start and end
#        cut -f1,2 {output.mp} | awk 'BEGIN{{
#        for(r=getline; r>0;){{
#            pos=$1;
#            for(s=e=$2; (r=getline)>0 && $2<=e+1; e=$2);
#            print s==e ? s : pos":"s-6"-"e+6
#        }}
#        exit -r
#        }}' > {output.bed)
#        cut -f1,2 {output.mp} | awk 'BEGIN{{
#        for(r=getline; r>0;){{
#            pos=$1;
#            for(s=e=$2; (r=getline)>0 && $2<=e+1; e=$2);
#            print s==e ? s : pos":"s-6"-"s-1"\n"pos":"e+1"-"e+6
#        }}
#        exit -r
#        }}' > {output.bed)
#        
#        samtools faidx /scratch/project_2001746/Racoon/references/Vulpes_vulpes.VulVul2.2.dna.toplevel.fa -r ./mockToRef.mpileup.bed -o flanking.fa