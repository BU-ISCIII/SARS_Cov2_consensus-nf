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

/*
 * SET UP CONFIGURATION VARIABLES
 */
params.help = false

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

/*
 * Default and custom value for configurable variables
 */

params.fasta = false
if( params.fasta ){
    fasta_file = file(params.fasta)
    if( !fasta_file.exists() ) exit 1, "Fasta file not found: ${params.fasta}."
}


// gtf file
params.gtf = false

if( params.gtf ){
    gtf_file = file(params.gtf)
    if( !gtf_file.exists() ) exit 1, "GTF file not found: ${params.gtf}."
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

if ( ! params.gtf ){
    exit 1, "GTF file not provided for assembly step, please declare it with --gtf /path/to/gtf_file"
}

/*
 * Create channel for input files
 */

// Create channel for input reads.
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { raw_reads_fastqc; raw_reads_trimming }

// Header log info
log.info "========================================="
log.info " BU-ISCIII/bacterial_wgs_training : WGS analysis practice v${version}"
log.info "========================================="
def summary = [:]
summary['Reads']               = params.reads
summary['Data Type']           = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Fasta Ref']           = params.fasta
summary['GTF File']            = params.gtf
summary['Keep Duplicates']     = params.keepduplicates
summary['Step']                = params.step
summary['Container']           = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']        = "$HOME"
summary['Current user']        = "$USER"
summary['Current path']        = "$PWD"
summary['Working dir']         = workflow.workDir
summary['Output dir']          = params.outdir
summary['Script dir']          = workflow.projectDir
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
