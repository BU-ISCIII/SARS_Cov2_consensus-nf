<img src="docs/images/BU_ISCIII_logo.png" alt="logo" width="300" align="right"/>

# SARS_Cov2_consensus-nf
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)

<!--
[![Docker](https://img.shields.io/docker/automated/nfcore/rnaseq.svg)](https://hub.docker.com/r/nfcore/rnaseq/)
-->
### Introduction

**BU-ISCIII/SARS_Cov2_consensus-nf** is a bioinformatics analysis pipeline used to analyze SARS-Cov2 Illumina data. The approach followed by this pipeline is to create a consensus genome through mapping, variant calling and consensus genome generation.

The workflow processes raw data from FastQ inputs ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)), maps the reads ([Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [Samtools](http://www.htslib.org/doc/samtools.html)) against the host and the viral reference genome. Optionally, if the data has been obtaine through amplicon sequencing you can use [iVar](http://gensoft.pasteur.fr/docs/ivar/1.0/manualpage.html) to trim the primers. Then the pipeline calls for variants and annotates them ([VarScan](http://varscan.sourceforge.net/), [SnpEff](http://snpeff.sourceforge.net/)), generates a genome sequence consensus ([Bgzip](http://www.htslib.org/doc/bgzip.html) and [BCFtools](http://www.htslib.org/doc/bcftools.html)) with the variants. Finally we generate a stats report with [MultiQC](https://multiqc.info/) See the [output documentation](docs/output.md) for more details of the results.

The pipeline is built using [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.

### Documentation
The BU-ISCIII/SARS_Cov2 pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. [Running the pipeline](docs/usage.md)
3. [Output and how to interpret the results](docs/output.md)
