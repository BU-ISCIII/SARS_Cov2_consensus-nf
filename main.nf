#!/usr/bin/env nextflow

/*
========================================================================================
                  SARS-Cov2 Illumina Consensus Analysis
========================================================================================
 #### Homepage / Documentation
 https://github.com/BU-ISCIII/SARS_Cov2_consensus-nf
 @#### Authors
 Sarai Varona <s.varona@isciii.es>
 Sara Monzón  <smonzon@isciii.es>
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Pipeline overview:
 - 1. : Preprocessing
 	- 1.1: FastQC - for raw sequencing reads quality control
 	- 1.2: Trimmomatic - raw sequence trimming
 - 2. : Mapping
 	- 2.1 : Bowtie2 - Mapping to host and reference viral genome
 	- 2.2 : Samtools - SAM and BAM files processing and stats
  - 2.3 : Picard - Mapping stats
  - 2.4 : iVar - Remove primers by coordinates.
 - 3. : Variant calling, annotation and consensus:
 	- 3.1 : VarScan - Variant calling
  - 3.2 : SnpEff - Variant Annotation
  - 3.3 : Bgzip - Variant calling vcf file compression
  - 3.4 : Bcftools - Consensus genome
 - 4. : Stats
  - 4.1 : MultiQC - Final html report with stats
 ----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
     BU-ISCIII/SARS_Cov2_consensus-nf : SARS_Cov2 Illumina consensus genome analysis v${version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run SARS_Cov2-nf/main.nf --reads '*_R{1,2}.fastq.gz' --viral_fasta ../../REFERENCES/NC_045512.2.fasta --viral_index '../REFERENCES/NC_045512.2.fasta.*' --host_fasta /processing_Data/bioinformatics/references/eukaria/homo_sapiens/hg38/UCSC/genome/hg38.fullAnalysisSet.fa --host_index '/processing_Data/bioinformatics/references/eukaria/homo_sapiens/hg38/UCSC/genome/hg38.fullAnalysisSet.fa.*' --amplicons_file ../REFERENCES/nCoV-2019.schemeMod.bed --outdir ./ -profile hpc_isciii

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes).
      --viral_fasta                 Path to Fasta reference
      --viral_index                 Path to viral fasta index
      --host_fasta                  Path to host Fasta sequence
      --host_index                  Path to host fasta index

    Options:
      --singleEnd                   Specifies that the input is single end reads
      --amplicons_file              Path to amplicons BED file.

    Trimming options
      --notrim                      Specifying --notrim will skip the adapter trimming step.
      --saveTrimmed                 Save the trimmed Fastq files in the the Results directory.
      --trimmomatic_adapters_file   Adapters index for adapter removal
      --trimmomatic_adapters_parameters Trimming parameters for adapters. <seed mismatches>:<palindrome clip threshold>:<simple clip threshold>. Default 2:30:10
      --trimmomatic_window_length   Window size. Default 4
      --trimmomatic_window_value    Window average quality requiered. Default 20
      --trimmomatic_mininum_length  Minimum length of reads

    Other options:
      --save_unmapped_host          Save the reads that didn't map to host genome
      --outdir                      The output directory where the results will be saved
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */
params.help = false


// Pipeline version
version = '1.0'

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

/*
 * Default and custom value for configurable variables
 */

params.viral_fasta = false
if( params.viral_fasta ){
    viral_fasta_file = file(params.viral_fasta)
    if( !viral_fasta_file.exists() ) exit 1, "Fasta file not found: ${params.viral_fasta}."
}

params.amplicons_file = false
if( params.amplicons_file ){
    amplicons_bed_file = file(params.amplicons_file)
    if ( !amplicons_bed_file.exists() ) exit 1, "Amplicons BAM file not found: $params.amplicons_file"
}

params.host_fasta = false
if( params.host_fasta ){
    host_fasta_file = file(params.host_fasta)
    if( !host_fasta_file.exists() ) exit 1, "Fasta file not found: ${params.host_fasta}."
}


// Output md template location
output_docs = file("$baseDir/docs/output.md")

// Trimming
// Trimming default
params.notrim = false
// Output files options
params.saveTrimmed = false
// Default trimming options
params.trimmomatic_adapters_file = "\$TRIMMOMATIC_PATH/adapters/NexteraPE-PE.fa"
params.trimmomatic_adapters_parameters = "2:30:10"
params.trimmomatic_window_length = "4"
params.trimmomatic_window_value = "20"
params.trimmomatic_mininum_length = "50"

// SingleEnd option
params.singleEnd = false

// Validate  mandatory inputs
params.reads = false
if (! params.reads ) exit 1, "Missing reads: $params.reads. Specify path with --reads"

/*
 * Create channel for input files
 */

// Create channel for input reads.
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { raw_reads_fastqc; raw_reads_trimming }

// Create channel for reference index
if( params.host_index ){
    Channel
        .fromPath(params.host_index)
        .ifEmpty { exit 1, "Host fasta index not found: ${params.host_index}" }
        .into { host_index_files }
}

if( params.viral_index ){
    Channel
        .fromPath(params.viral_index)
        .ifEmpty { exit 1, "Viral fasta index not found: ${params.viral_index}" }
        .into { viral_index_files; viral_index_files_variant_calling }
}


// Header log info
log.info "========================================="
log.info " BU-ISCIII/bacterial_wgs_training : WGS analysis practice v${version}"
log.info "========================================="
def summary = [:]
summary['Reads']               = params.reads
summary['Data Type']           = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Fasta Ref']           = params.viral_fasta
summary['Container']           = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']        = "$HOME"
summary['Current user']        = "$USER"
summary['Current path']        = "$PWD"
summary['Working dir']         = workflow.workDir
summary['Output dir']          = params.outdir
summary['Script dir']          = workflow.projectDir
summary['Save Unmapped']        = params.save_unmapped_host
summary['Save Trimmed']        = params.saveTrimmed
if( params.notrim ){
    summary['Trimming Step'] = 'Skipped'
} else {
    summary['Trimmomatic adapters file'] = params.trimmomatic_adapters_file
    summary['Trimmomatic adapters parameters'] = params.trimmomatic_adapters_parameters
    summary["Trimmomatic window length"] = params.trimmomatic_window_length
    summary["Trimmomatic window value"] = params.trimmomatic_window_value
    summary["Trimmomatic minimum length"] = params.trimmomatic_mininum_length
}
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "===================================="

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '0.25.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

/*
 * STEP 1.1 - FastQC
 */
process fastqc {
	tag "$prefix"
	publishDir "${params.outdir}/01-fastQC", mode: 'copy',
		saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

	input:
	set val(name), file(reads) from raw_reads_fastqc

	output:
	file '*_fastqc.{zip,html}' into fastqc_results
	file '.command.out' into fastqc_stdout

	script:

	prefix = name - ~/(_S[0-9]{2})?(_L00[1-9])?(.R1)?(_1)?(_R1)?(_trimmed)?(_val_1)?(_00*)?(\.fq)?(\.fastq)?(\.gz)?$/
	"""
	fastqc -t 1 $reads
	"""
}

/*
 * STEPS 1.2 Trimming
 */
process trimming {
	tag "$prefix"
	publishDir "${params.outdir}/02-preprocessing", mode: 'copy',
		saveAs: {filename ->
			if (filename.indexOf("_fastqc") > 0) "../03-preprocQC/$filename"
			else if (filename.indexOf(".log") > 0) "logs/$filename"
      else if (params.saveTrimmed && filename.indexOf(".fastq.gz")) "trimmed/$filename"
			else null
	}

	input:
	set val(name), file(reads) from raw_reads_trimming

	output:
	file '*_paired_*.fastq.gz' into trimmed_paired_reads,trimmed_paired_reads_bwa,trimmed_paired_reads_bwa_virus
	file '*_unpaired_*.fastq.gz' into trimmed_unpaired_reads
	file '*_fastqc.{zip,html}' into trimmomatic_fastqc_reports
	file '*.log' into trimmomatic_results

	script:
	prefix = name - ~/(_S[0-9]{2})?(_L00[1-9])?(.R1)?(_1)?(_R1)?(_trimmed)?(_val_1)?(_00*)?(\.fq)?(\.fastq)?(\.gz)?$/
	"""
	java -jar $TRIMMOMATIC_PATH/trimmomatic-0.33.jar PE -threads 1 -phred33 $reads $prefix"_paired_R1.fastq" $prefix"_unpaired_R1.fastq" $prefix"_paired_R2.fastq" $prefix"_unpaired_R2.fastq" ILLUMINACLIP:${params.trimmomatic_adapters_file}:${params.trimmomatic_adapters_parameters} SLIDINGWINDOW:${params.trimmomatic_window_length}:${params.trimmomatic_window_value} MINLEN:${params.trimmomatic_mininum_length} 2> ${name}.log

	gzip *.fastq

	fastqc -q *_paired_*.fastq.gz

	"""
}
/*
 * STEPS 2.1 Mapping host
 */
process mapping_host {
	tag "$prefix"
	publishDir "${params.outdir}/04-mapping_host", mode: 'copy',
		saveAs: {filename ->
			if (filename.indexOf(".bam") > 0) "mapping/$filename"
			else if (filename.indexOf(".bai") > 0) "mapping/$filename"
      else if (filename.indexOf(".txt") > 0) "stats/$filename"
      else if (filename.indexOf(".stats") > 0) "stats/$filename"
	}

	input:
	set file(readsR1),file(readsR2) from trimmed_paired_reads_bwa
  file refhost from host_fasta_file
  file index from host_index_files.collect()

	output:
	file '*_sorted.bam' into mapping_host_sorted_bam
  file '*.bam.bai' into mapping_host_bai
	file '*_flagstat.txt' into mapping_host_flagstat
	file '*.stats' into mapping_host_picardstats

	script:
	prefix = readsR1.toString() - '_paired_R1.fastq.gz'
	"""
  bowtie2 --local -x $refhost -1 $readsR1 -2 $readsR2 --very-sensitive-local -S $prefix".sam"
  samtools sort -o $prefix"_sorted.bam" -O bam -T $prefix $prefix".sam"
  samtools index $prefix"_sorted.bam"
  samtools flagstat $prefix"_sorted.bam" > $prefix"_flagstat.txt"
  picard CollectWgsMetrics COVERAGE_CAP=1000000 I=$prefix"_sorted.bam" O=$prefix".stats" R=$refhost
	"""
}

/*
 * STEPS 2.2 Mapping virus
 */
process mapping_virus {
	tag "$prefix"
	publishDir "${params.outdir}/05-mapping_virus", mode: 'copy',
		saveAs: {filename ->
			if (filename.indexOf(".bam") > 0) "mapping/$filename"
			else if (filename.indexOf(".bai") > 0) "mapping/$filename"
      else if (filename.indexOf(".txt") > 0) "stats/$filename"
      else if (filename.indexOf(".stats") > 0) "stats/$filename"
	}

	input:
	set file(readsR1),file(readsR2) from trimmed_paired_reads_bwa_virus
  file refvirus from viral_fasta_file
  file index from viral_index_files.collect()

	output:
	file '*_sorted.bam' into mapping_virus_sorted_bam,mapping_virus_sorted_bam_variant_calling
  file '*.bam.bai' into mapping_virus_bai,mapping_virus_bai_variant_calling
	file '*_flagstat.txt' into mapping_virus_flagstat
	file '*.stats' into mapping_virus_picardstats

	script:
  prefix = readsR1.toString() - '_paired_R1.fastq.gz'
	"""
  bowtie2 --local -x $refvirus -1 $readsR1 -2 $readsR2 --very-sensitive-local -S $prefix".sam"
  samtools sort -o $prefix"_sorted.bam" -O bam -T $prefix $prefix".sam"
  samtools index $prefix"_sorted.bam"
  samtools flagstat $prefix"_sorted.bam" > $prefix"_flagstat.txt"
  picard CollectWgsMetrics COVERAGE_CAP=1000000 I=$prefix"_sorted.bam" O=$prefix".stats" R=$refvirus
	"""
}

/*
 * STEPS 2.3 iVar trim primers
 */
if (params.amplicons_file) {
  process ivar_trimming {
	 tag "$prefix"
	  publishDir path: { "${params.outdir}/05-mapping_virus/ivar" }, mode: 'copy'

	  input:
	  file sorted_bam from mapping_virus_sorted_bam_variant_calling
    file amplicons_bed from amplicons_bed_file

	  output:
	  file '*_primertrimmed_sorted.bam' into ivar_sorted_bam,sorted_bam_variant_calling,ivar_sorted_bam_consensus
    file '*_primertrimmed_sorted.bam.bai' into ivar_bai,bam_bai_variant_calling,ivar_bai_consensus

	  script:
    prefix = sorted_bam.baseName - ~/(_S[0-9]{2})?(_L00[1-9])?(.R1)?(_1)?(_R1)?(_sorted)?(_paired)?(_00*)?(\.bam)?(\.fastq)?(\.gz)?$/
	  """
    samtools view -b -F 4 $sorted_bam > $prefix"_onlymapped.bam"
    samtools index $prefix"_onlymapped.bam"
    ivar trim -e -i $prefix"_onlymapped.bam" -b $amplicons_bed -p $prefix"_primertrimmed" -q 15 -m 50 -s 4
    samtools sort -o $prefix"_primertrimmed_sorted.bam" -O bam -T $prefix $prefix"_primertrimmed.bam"
    samtools index $prefix"_primertrimmed_sorted.bam"
	  """
  }
} else {
    mapping_virus_sorted_bam_variant_calling
      .set {sorted_bam_variant_calling}
    mapping_virus_bai_variant_calling
      .set {bam_bai_variant_calling}
}


/*
 * STEPS 3.1 Variant Calling
 */
process variant_calling {
	tag "$prefix"
	publishDir "${params.outdir}/06-variant_calling", mode: 'copy',
		saveAs: {filename ->
			if (filename.indexOf(".pileup") > 0) "pileup/$filename"
			else if (filename.indexOf("_majority.vcf") > 0) "majority_allele/$filename"
      else if (filename.indexOf(".vcf") > 0) "lowfreq_vars/$filename"
			else params.saveTrimmed ? filename : null
	}

	input:
	file sorted_bam from sorted_bam_variant_calling
  file bam_index from bam_bai_variant_calling
  file refvirus from viral_fasta_file
  file index from viral_index_files_variant_calling.collect()

	output:
	file '*.pileup' into variant_calling_pileup
  file '*_majority.vcf' into majority_allele_vcf,majority_allele_vcf_annotation,majority_allele_vcf_consensus
	file '*_lowfreq.vcf' into lowfreq_variants_vcf,lowfreq_variants_vcf_annotation

	script:
	prefix = sorted_bam.baseName - ~/(_primertrimmed)?(_L00[1-9])?(.R1)?(_1)?(_R1)?(_sorted)?(_paired)?(_00*)?(\.bam)?(\.fastq)?(\.gz)?$/
	"""
  samtools mpileup -A -d 20000 -Q 0 -f $refvirus $sorted_bam > $prefix".pileup"
  varscan mpileup2cns $prefix".pileup" --min-var-freq 0.02 --p-value 0.99 --variants --output-vcf 1 > $prefix"_lowfreq.vcf"
  varscan mpileup2cns $prefix".pileup" --min-var-freq 0.8 --p-value 0.05 --variants --output-vcf 1 > $prefix"_majority.vcf"
	"""
}

/*
 * STEPS 3.2 Variant Calling annotation
 */
process variant_calling_annotation {
 	tag "$prefix"
  publishDir "${params.outdir}/07-annotation", mode: 'copy',
		saveAs: {filename ->
			if (filename.indexOf("majority.ann.vcf") > 0) "majority/$filename"
			else if (filename.indexOf("majority_snpEff_genes.txt") > 0) "majority/$filename"
      else if (filename.indexOf("majority_snpEff_summary.html") > 0) "majority/$filename"
      else if (filename.indexOf("lowfreq.ann.vcf") > 0) "lowfreq/$filename"
      else if (filename.indexOf("lowfreq_snpEff_genes.txt") > 0) "lowfreq/$filename"
      else if (filename.indexOf("lowfreq_snpEff_summary.html") > 0) "lowfreq/$filename"
	}

 	input:
 	file majority_variants from majority_allele_vcf_annotation
  file low_variants from lowfreq_variants_vcf_annotation

 	output:
 	file '*_majority.ann.vcf' into majority_annotated_variants
  file '*_majority_snpEff_genes.txt' into majority_snpeff_genes
 	file '*_majority_snpEff_summary.html' into majority_snpeff_summary
  file '*_lowfreq.ann.vcf' into lowfreq_annotated_variants
  file '*_lowfreq_snpEff_genes.txt' into lowfreq_snpeff_genes
 	file '*_lowfreq_snpEff_summary.html' into lowfreq_snpeff_summary

 	script:
 	prefix = majority_variants.baseName - ~/(_S[0-9]{2})?(_majority)?(.R1)?(_1)?(_R1)?(_sorted)?(_paired)?(_00*)?(\.bam)?(\.vcf)?(\.gz)?$/
 	"""
  snpEff sars-cov-2 $majority_variants > $prefix"_majority.ann.vcf"
  mv snpEff_genes.txt $prefix"_majority_snpEff_genes.txt"
  mv snpEff_summary.html $prefix"_majority_snpEff_summary.html"
  snpEff sars-cov-2 $low_variants > $prefix"_lowfreq.ann.vcf"
  mv snpEff_genes.txt $prefix"_lowfreq_snpEff_genes.txt"
  mv snpEff_summary.html $prefix"_lowfreq_snpEff_summary.html"
 	"""
}

/*
 * STEPS 3.3 Consensus Genome
 */
process genome_consensus {
  tag "$prefix"
  publishDir "${params.outdir}/08-mapping_consensus", mode: 'copy',
		saveAs: {filename ->
			if (filename.indexOf("_consensus.fasta") > 0) "consensus/$filename"
			else if (filename.indexOf("_consensus_masked.fasta") > 0) "masked/$filename"
	}

  input:
  file variants from majority_allele_vcf_consensus
  file refvirus from viral_fasta_file
  file ivar_bam from ivar_sorted_bam_consensus
  file ivar_bai from ivar_bai_consensus

  output:
  file '*_consensus.fasta' into consensus_fasta
  file '*_consensus_masked.fasta' into masked_fasta

  script:
  prefix = variants.baseName - ~/(_majority)?(_paired)?(\.vcf)?(\.gz)?$/
  refname = refvirus.baseName - ~/(\.2)?(\.fasta)?$/
  """
  bgzip -c $variants > $prefix"_"$refname".vcf.gz"
  bcftools index $prefix"_"$refname".vcf.gz"
  cat $refvirus | bcftools consensus $prefix"_"$refname".vcf.gz" > $prefix"_"$refname"_consensus.fasta"
  bedtools genomecov -bga -ibam $ivar_bam -g $refvirus | awk '\$4 < 20' | bedtools merge > $prefix"_"$refname"_bed4mask.bed"
  bedtools maskfasta -fi $prefix"_"$refname"_consensus.fasta" -bed $prefix"_"$refname"_bed4mask.bed" -fo $prefix"_"$refname"_consensus_masked.fasta"
  sed -i 's/NC_045512.2/$prefix/g' $prefix"_"$refname"_consensus_masked.fasta"
  """
}
