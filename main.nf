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
