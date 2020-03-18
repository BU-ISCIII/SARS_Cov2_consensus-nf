<img src="docs/images/BU_ISCIII_logo.png" alt="logo" width="300" align="right"/>

# SARS_Cov2-nf
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)

<!--
[![Docker](https://img.shields.io/docker/automated/nfcore/rnaseq.svg)](https://hub.docker.com/r/nfcore/rnaseq/)
-->
### Introduction

**BU-ISCIII/SARS_Cov2-nf** is a bioinformatics analysis pipeline used to analyze SARS-Cov2 Illumina SISPA data.

The workflow processes raw data from FastQ inputs ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)), maps the reads ([BWA]() and [Samtools]()) against the host and the viral reference genome, calls for variants and annotates them ([VarScan](), [SnpEff]()), generates a genome sequence consensus ([Bgzip]() and [BCFtools]()) with the variants, generated a _De Novo_ genome assembly ([Spades]() and/or [Unicycler]()), reorders the contigs generating a draft genome ([Abacas]()), aligns the assemblies ([BLAST]()) and creates a series of stats and graphs ([plasmidID]()). See the [output documentation](docs/output.md) for more details of the results.

The pipeline is built using [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.

<!--
### Documentation
The BU-ISCIII/SARS_Cov2 pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
    * [Reference genomes](docs/configuration/reference_genomes.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)
-->
