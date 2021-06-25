# Snakemake workflow: snakemake_rnaseq_combined_workflow

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.2.1-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/rna-seq-star-deseq2.svg?branch=master)](https://travis-ci.org/snakemake-workflows/rna-seq-star-deseq2)
[![Snakemake-Report](https://img.shields.io/badge/snakemake-report-green.svg)](https://cdn.rawgit.com/snakemake-workflows/rna-seq-star-deseq2/master/.test/report.html)

This workflow performs trimming with trimmomatic sliding window, runs alignment with STAR and/or HISAT2, produces a featureCounts file and then runs QC before and after trimming to generate 3 multiqc reports..

## Authors

* Hans Vasquez-Gross (@hansvg)

## Usage

### Simple

#### Step 1: Install workflow

clone this workflow to your local computer

#### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml` and units.tsv.

To help fill out the units.tsv, you can use MS Excel to create the 4 columns sample,units,fq1,fq2

Then in the directory where you have your fastq files, run the following command to get the full path to the files. This can be easily copy/pasted into the correct column.
`ls -d1 $PWD/*R1.fastq.gz`

Likewise, once the worksheet is filled out, copy the table and paste it into an empty units.tsv file.

#### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

