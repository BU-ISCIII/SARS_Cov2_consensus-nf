# Output description for SARS-Cov2 consensus genome pipeline
**SARS_Cov2_consensus-nf** is a bioinformatics best-practice analysis pipeline used for SARS-Cov2 consensus genome generation. The pipeline is focused in mapping and variants calling.

**Note**:If you want this PDF file to work you can't move it from DOC folder and you can't either move the files inside ANALYSIS folder.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* <a href="#fastqc">FastQC</a> v0.11.8 - read quality control
* <a href="#trimming">Trimmomatic</a> v.0.38 - adapter and low quality trimming.
* <a href="#bowtie">Bowtie</a> v2.3.5 - mapping against reference genomes.
* <a href="#samtools">SAMtools</a> v1.2 - Mapping result processing and unmapped reads selection.
* <a href="#picard">Picard</a> v1.140 - Enrichment and alignment metrics.
* <a href="#ivar">iVar</a> v1.1 - Amplicons primers trimming by position.
* <a href="#varscan">VarScan</a> v2.4.4-0 - Variant calling.
* <a href="#snpeff">snpEff and SnpSift</a> 4.3.1t - Variant calling annotation.
* <a href="#bcftools">Bcftools</a> v1.9 - Variant calling index and consensus genome generation.
* <a href="#bedtools">Bedtools</a> v2.26.0 - Consensus genome masking.
* <a href="#multiqc">MultiQC</a> v1.7 - aggregate report, describing results of the whole pipeline.

**Reference genome:** NC_045512.2

Depending on the analysis, we will have some ANALYSIS_IDs. This ANALYSIS_IDs are going to be composed of the date of the analysis, and some identification of the type of analysis.

## Preprocessing
### FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)<a name="fastqc"> </a><a href="#fastqc_reference">[1]</a> gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output directory: `01-fastQC`**

* `{sample_id}_R[12]_fastqc.html`
  * html report. This file can be opened in your favourite web browser (Firefox/chrome preferable) and it contains the different graphs that fastqc calculates for QC.
* `zips/{sample_id}_R[12]_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

### Trimming
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)<a name="trimming"> </a><a href="#trimming_reference">[2]</a> is used for removal of adapter contamination and trimming of low quality regions.
Parameters included for trimming are:
-  Nucleotides with phred quality < 10 in 3'end.
-  Mean phred quality < 20 in a 4 nucleotide window.
-  Read lenght < 50

**Results directory: `02-preprocessing`**
- Files:
   - `trimmed/{sample_id}_R[12]_filtered.fastq.gz`: contains high quality reads with both forward and reverse tags surviving.
   - `trimmed/{sample_id}_R[12]_unpaired.fastq.gz`: contains high quality reads with only forward or reverse tags surviving.
   - `logs/{sample_id}.log`: trimming log file.

 **Note**:To see how your reads look after trimming, look at the FastQC reports in the ANALYSIS/{ANALYSIS_ID}/03-preprocQC directory

 **Note**:From now on, all the steps will be host specific.

## Mapping and Primer trimming
### Bowtie and Samtools
[Bowtie](http://bio-bwa.sourceforge.net/)<a name="bowtie"> </a><a href="#Bowtie_reference">[3]</a> is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of about 50 up to 100s of characters to relatively long genomes. Bowtie 2 indexes the genome with an FM Index (based on the Burrows-Wheeler Transform or BWT) to keep its memory footprint small. Bowtie 2 supports gapped, local, and paired-end alignment modes.

The result mapping files are further processed with [SAMtools](http://samtools.sourceforge.net/)<a name="samtools"> </a><a href="#SAMtools_reference">[4]</a>, sam format is converted to bam, sorted and an index .bai is generated.

We mapped the fastq file agains both reference host genome and reference viral genome.

**Output directory: `0[4/5]-mapping_[host/virus]`**

* `mapping/{sample_id}_sorted.bam`
  * Sorted aligned bam file.
* `mapping/{sample_id}_sorted.bam.bai`
  * Index file for soreted aligned bam.
* `stats/{sample_id}_flagstat.txt`
  * Samtools flagstats mapping stats summary.

### Picard
[Picard](https://broadinstitute.github.io/picard/index.html)<a name="picard"> </a><a href="#Picard_reference">[5]</a> is a set of command line tools for manipulating high-throughput sequencing (HTS) data. In this case we used it to obtain mapping stats.

**Output directory: `0[4/5]-mapping_[host/virus]`**

* `stats/{sample_id}.stats`
  * Picard metrics summary file for evaluating coverage and performance.

Picard documentation: [Picarddocs](https://broadinstitute.github.io/picard/command-line-overview.html)

### iVar
[iVar](http://gensoft.pasteur.fr/docs/ivar/1.0/manualpage.html)<a name="ivar"> </a><a href="#iVar_reference">[6]</a> uses primer positions supplied in a BED file to soft clip primer sequences from an aligned, sorted and indexed BAM file. Following this, the reads are trimmed based on a quality threshold through a sliding window approach that slides from the 5' end to the 3' end and if at any point the average base quality in the window falls below the threshold, the remaining read is soft clipped. If after trimming, the length of the read is greater than the minimum length specified, the read is written to the new trimmed BAM file.

**Output directory: `05-mapping_virus`**

* `{sample_id}_primertrimmed_sorted.bam`
  * Sorted aligned bam file after trimming.
* `{sample_id}_primertrimmed_sorted.bam.bai`
  * Index file for sorted aligned trimmed bam.
* `{sample_id}_primertrimmed_flagstats.txt`
  * Samtools flagstats summary file.
* `{sample_id}_primertrimmed.stats`
  * Picard metrics summary file for evaluating coverage and performance.

## Variant calling
In this pipeline to generate the consensus viral genome we use the approach of calling for variants between the mapped reads and the reference viral genome, and adding these variants to the reference viral genome.

### Variant calling
#### VarScan
First of all SAMtools is used to generate the variant calling VCF file. Then [VarScan](http://varscan.sourceforge.net/)<a name="varscan"> </a><a href="#VarScan_reference">[7]</a> is used to call for major and low frequency variants. VarScan is a platform-independent software tool developed at the Genome Institute at Washington University to detect variants in NGS data.

**Output directory: `06-variant_calling`**

* `pileup/{sample_id}.pileup`
  * Variant calling pileup file. The pileup file summarises all data from the reads at each genomic region that is covered by at least one read. Each row of the pileup file gives similar information to a single vertical column of reads in the IGV view.
* `lowfreq_vars/{sample_id}_lowfreq.vcf`
  * Low frequency variants VCF.
* `majority_allele/{sample_id}_majority.vcf`
  * Majority variants VCF file.

### Annotation
#### SnpEff and SnpSift
[SnpEff](http://snpeff.sourceforge.net/SnpEff.html)<a name="snpeff"> </a><a href="#SnpEff_reference">[8]</a> is a genetic variant annotation and functional effect prediction toolbox. It annotates and predicts the effects of genetic variants on genes and proteins (such as amino acid changes). [SnpSift](http://snpeff.sourceforge.net/SnpSift.html) annotates genomic variants using databases, filters, and manipulates genomic annotated variants. Once you annotated your files using SnpEff, you can use SnpSift to help you filter large genomic datasets in order to find the most significant variants for your experiment.

**Output directory: `07-annotation`**

* `lowfreq/{sample_id}_lowfreq.ann.table.txt`
  * Low frequency variants SnpSift summary table.
* `lowfreq/{sample_id}_lowfreq.ann.vcf`
  * Low frequency variants annotated VCF table.
* `lowfreq/{sample_id}_lowfreq_snpEff_genes.txt`
  * Low frequency variants genes table.
* `lowfreq/{sample_id}_lowfreq_snpEff_summary.html`
  * Low frequency variants summary html file.
* `majority/{sample_id}_majority.ann.table.txt`
  * Majority variants SnpSift summary table.
* `majority/{sample_id}_majority.ann.vcf`
  * Majority variants annotated VCF table.
* `majority/{sample_id}_majority_snpEff_genes.txt`
  * Majority variants genes table.
* `majority/{sample_id}_majority_snpEff_summary.html`
  * Majority variants summary html file.

## Consensus genome
### BCFtools
[Bcftools](http://samtools.github.io/bcftools/bcftools.html)<a name="bcftools"> </a><a href="#BCFtools_reference">[9]</a> is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. The resulting variant calling vcf for haploid genomes is indexed and then the consensus genome is created adding the variants to the reference viral genome. This consensus genome was obtained using the predominant variants (majority) of the mapping file.

**Output directory: `08-mapping_consensus`**

* `consensus/{sample_id}_{reference_virus_name}_consensus.fasta`
  * Consensus viral genome file generated from adding the variants called before to the viral reference genome. These variants are only the majoritarian variants, inlcuding only SNPs and small indels. This file is also contained in the 10-final_results folder as {sample_id}_{reference_virus_name}_consensus.fasta.

### Bedtools
[Bedtools](https://bedtools.readthedocs.io/en/latest/)<a name="bedtools"> </a><a href="#Bedtools_reference">[10]</a> are a swiss-army knife of tools for a wide-range of genomics analysis tasks. In this case we use:
  * bedtools genomecov computes histograms (default), per-base reports (-d) and BEDGRAPH (-bg) summaries of feature coverage (e.g., aligned sequences) for a given genome.
  * bedtools maskfasta masks sequences in a FASTA file based on intervals defined in a feature file. This may be useful fro creating your own masked genome file based on custom annotations or for masking all but your target regions when aligning sequence data from a targeted capture experiment.

**Output directory: `08-mapping_consensus`**

* `masked/{sample_id}_{reference_virus_name}_consensus_masked.fasta`
  * Masked consensus fasta file.

## Stats
### MultiQC
[MultiQC](http://multiqc.info)<a name="multiqc"> </a><a href="#MultiQC_reference">[11]</a> is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualized in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `99-stats/MultiQC/`**

* `multiqc_report.html`
  * MultiQC report: standalone HTML file that can be viewed in your web browser
* `multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)

## Final Results
We hace collected the most significant files for you.

**Output directory:** `RESULTS`

* draft_genomes: this folder contains the draft genomes.
* ordered_contigs: this folder contains the ordered contigs for each sample.
* circos_images: circos images for the reconstructed genomes.
* reads_stats: statistics for mapped reads against host and virus, with coments.
