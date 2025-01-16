# High-Power Computer Genotyping-by-Sequencing Pipeline

## Introduction

This git repository houses scripts that can be utilized on [SCINet](https://scinet.usda.gov/) computing clusters to perform the [TASSEL 5 Standalone](https://tassel.bitbucket.io/) genotyping-by-sequencing (GBS) pipeline. This pipeline is inteded for educational use only and comes with ***NO WARRENTY***.

## Overall Workflow

This TASSEL pipeline has four distinct sections:

### Discovery

The discovery step involves making a [sqlite](https://www.sqlite.org/index.html) database of all possible single nucleotide polymorphims (SNPs) of FASTQ files aligned to a reference genome of the user's choice.

### Production

The production step involves making a variant call format ([VCF](https://en.wikipedia.org/wiki/Variant_Call_Format)) file for all individuals in a requested panel. 

### Filtration 

The filtration step involves using parameters like read depth and missing data to filter the raw production VCF.

### Imputation

The imputation step involves using the [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html) algorithem for imputation of the filtered VCF.

## Instructions

The pipeline is currently set up to be run on the [Atlas HPC](https://www.hpc.msstate.edu/computing/atlas/) which is a collaborator with SCINet. This script takes advantage of the [environment-module system](https://modules.sourceforge.net/) which allows for the selective loading of required packages. Atlas' module of TASSEL appears to be defunct at the moment, therefore this script pull the TASSEL 5 Standalone package directly from GitHub for use in this pipeline. 

## Modifying the Header of the .slurm Script

The main script used in this pipeline is `atlas_gbs_pipeline.slurm`
