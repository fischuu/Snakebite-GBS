library(DiagrammeR)

# Create a DiagrammeR graph
graph <- grViz("
  digraph snakemake_pipeline {
    graph [layout = dot, rankdir = LR]

    node [shape = box, style = filled, fillcolor = lightblue]
    edge [arrowhead = vee]

    # Define new reads nodes
      input_reads [label = 'Input Reads (fastq)', shape = folder]
      ref_genome [label = 'Reference genome (fasta)', shape = folder]
      insilico_genome [label = 'In-silico genome (fasta)', shape = folder]
      mockref_genome [label = 'Mock reference (fasta)', shape = folder]
      finalmock_genome [label = 'Final mock (fasta)', shape = folder]
    

    # Preparation
      concatenate_reads [label = 'Concatenate reads \n cat']
      trim_reads [label = 'Trim reads \n trimmomatic']
      insilico_digestion [label = 'In-silico digestion \n R']

    # Mock ref
      merge_reads [label = 'Merge reads \n PEAR']
      stitch_reads [label = 'Stich unmergable reads \n perls']
      assemble_genome [label = 'Assemble genome \n vsearch']
      refine_mock [label = 'Refine mock \n R']
      
    # Alignment
      align_reads_to_insilico [label = 'Align reads to in-silico \n BWA_mem']
      align_reads_to_ref [label = 'Align reads to ref \n BWA_mem']
      align_reads_to_mock [label = 'Align reads to mock \n BWA_mem']
      align_reads_to_final [label = 'Align reads to final mock \n BWA_mem']
      
    # Call variants
      create_pileup_ref [label = 'Create MPILEUP (Reference) \n samtools']
      create_pileup_mock [label = 'Create MPILEUP (Mock ref) \n samtools']
      create_pileup_final [label = 'Create MPILEUP (Final mock) \n samtools']
      create_mm_ref [label = 'Create Mastermatric (Reference) \n perl']
      create_mm_mock [label = 'Create Mastermatric (Mock ref) \n perl']
      create_mm_final [label = 'Create Mastermatric (Final mock) \n perl']
      filter_variants_ref [label = 'Filter variants (Reference) \n perl']
      filter_variants_mock [label = 'Filter variants (Mock ref) \n perl']
      filter_variants_final [label = 'Filter variants (Final ref) \n perl']
      create_vcf_ref [label = 'Create vcf (Reference) \n perl']
      create_vcf_mock [label = 'Create vcf (Mock ref) \n perl']
      create_vcf_final [label = 'Create vcf (Final mock) \n perl']
      
    # Reporting
      create_report [label = 'Create report \n R']

################################################################################
##
## Edges
##

    # Define edges for input reads
      input_reads -> concatenate_reads
      concatenate_reads -> trim_reads

    # Insilico
      ref_genome -> insilico_digestion
      insilico_digestion -> insilico_genome
      
    # Mock ref
      trim_reads -> merge_reads
      trim_reads -> stitch_reads
      merge_reads -> assemble_genome
      stitch_reads -> assemble_genome
      assemble_genome -> mockref_genome      
      
    # Alignments
      trim_reads -> align_reads_to_ref
      ref_genome -> align_reads_to_ref
      trim_reads -> align_reads_to_mock
      mockref_genome -> align_reads_to_mock
      trim_reads -> align_reads_to_final
      finalmock_genome -> align_reads_to_final
      trim_reads -> align_reads_to_insilico
      insilico_genome -> align_reads_to_insilico
      mockref_genome -> refine_mock
      align_reads_to_mock -> refine_mock
      refine_mock -> finalmock_genome
      
  # Variant calling
     align_reads_to_ref -> create_pileup_ref
     create_pileup_ref -> create_mm_ref
     create_mm_ref -> filter_variants_ref
     filter_variants_ref -> create_vcf_ref
     
     align_reads_to_mock -> create_pileup_mock
     create_pileup_mock -> create_mm_mock
     create_mm_mock -> filter_variants_mock
     filter_variants_mock -> create_vcf_mock
     
     align_reads_to_final -> create_pileup_final
     create_pileup_final -> create_mm_final
     create_mm_final -> filter_variants_final
     filter_variants_final -> create_vcf_final
     
     create_vcf_ref -> create_report
     create_vcf_mock -> create_report
     create_vcf_final -> create_report
     
################################################################################
##
## subgraphs
##

    subgraph cluster_preprocess_reads {
      label = 'Preprocessing'
      style = dashed
      color = lightgray

      concatenate_reads
      trim_reads
    }
    
    subgraph cluster_mockref {
      label = 'Mock reference'
      style = dashed
      color = lightgray

      merge_reads
      stitch_reads
      assemble_genome
    }
    
    subgraph cluster_refinemock {
      label = 'Refine Mock'
      style = dashed
      color = lightgray

      refine_mock
    }
    
    subgraph cluster_alignments {
      label = 'Alignment'
      style = dashed
      color = lightgray

      align_reads_to_ref
      align_reads_to_mock
      align_reads_to_final
    }
    
    subgraph cluster_insilico {
      label = 'In-silico digestion'
      style = dashed
      color = lightgray

      insilico_digestion
      align_reads_to_insilico
    }
  
    subgraph cluster_genomes {
      label = 'Genomes'
      style = dashed
      color = lightgray

      ref_genome
      mockref_genome
      finalmock_genome
      insilico_genome
    }
      
    
    subgraph cluster_variant_calling {
      label = 'Variant calling'
      style = dashed
      color = lightgray

      create_pileup_ref
      create_pileup_mock
      create_pileup_final
      create_mm_ref
      create_mm_mock
      create_mm_final
      filter_variants_ref
      filter_variants_mock
      filter_variants_final
      create_vcf_ref
      create_vcf_mock
      create_vcf_final
    }
    
  }
")

# Render the graph
graph

# Export to SVG
# svg_code <- export_svg(graph)

# Save as an HTML file
#html_code <- paste0('<html><body>', svg_code, '</body></html>')
#write(html_code, file = "pipeline_diagram.html")
