#!/usr/bin/env nextflow

/*
========================================================================================
                  SARS-Cov2 Illumina SISPA Analysis
========================================================================================
 #### Homepage / Documentation
 https://github.com/BU-ISCIII/SARS_Cov2-nf
 @#### Authors
 Sarai Varona <s.varona@isciii.es>
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Pipeline overview:
 - 1. : Preprocessing
 	- 1.1: FastQC - for raw sequencing reads quality control
 	- 1.2: Trimmomatic - raw sequence trimming
 - 2. : Mapping
 	- 2.1 : BWA - Mapping to host and reference viral genome
 	- 2.2 : Samtools - SAM and BAM files processing and stats
  - 2.3 : Picard - Mapping stats
 - 3. : Variant calling, annotation and consensus:
 	- 3.1 : VarScan - Variant calling
  - 3.2 : SnpEff - Variant Annotation
  - 3.3 : Bgzip - Variant calling vcf file compression
  - 3.4 : Bcftools - Consensus genome
 - 4. : DeNovo assembly:
  - 4.1 : Spades - De Novo assembly (normal mode and metaSpades mode)
  - 4.2 : Unicycler - De novo assembly
  - 4.3 : Quast - Assembly quality assessment.
  - 4.4 : ABACAS - Assembly contig reordering and draft generation
 - 5. : Assembly Alignment
  - 5.1 : BLAST - Assembly alignment to viral reference
 - 6. Stats & Graphs :
  - 6.1 : PlasmidID - Assembly plot generation
 ----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
     BU-ISCIII/SARS_Cov2-nf : SARS_Cov2 Illumina SISPA data analysis v${version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run SARS_Cov2-nf/main.nf --reads '*_R{1,2}.fastq.gz' --fasta ../../REFERENCES/NC_045512.2.fasta --gff ../../REFERENCES/NC_045512.2.gff

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes).
      --fasta                       Path to Fasta reference
      --gff					            		Path to GFF reference file. (Mandatory if step = assembly)

    Options:
      --singleEnd                   Specifies that the input is single end reads

    Trimming options
      --notrim                      Specifying --notrim will skip the adapter trimming step.
      --saveTrimmed                 Save the trimmed Fastq files in the the Results directory.
      --trimmomatic_adapters_file   Adapters index for adapter removal
      --trimmomatic_adapters_parameters Trimming parameters for adapters. <seed mismatches>:<palindrome clip threshold>:<simple clip threshold>. Default 2:30:10
      --trimmomatic_window_length   Window size. Default 4
      --trimmomatic_window_value    Window average quality requiered. Default 20
      --trimmomatic_mininum_length  Minimum length of reads

    Other options:
      --outdir                      The output directory where the results will be saved
    """.stripIndent()
}
