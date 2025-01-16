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

The main script used in this pipeline is [atlas_tassel_pipeline.slurm](https://github.com/zjwinn/HPC-GBS-Pipeline/blob/main/atlas_tassel_pipeline.slurm). The header of this script is formatted to the system specifications of atlas and looks like the following

```bash
#!/bin/bash
#SBATCH --job-name="tassel_pl"                   
#SBATCH --qos=normal
#SBATCH -p atlas                                 
#SBATCH -A guedira_seq_map                       
#SBATCH -N 1                                      
#SBATCH -n 48                                     
#SBATCH -t 7-00:00:00                             
#SBATCH --mail-user=YOUR.EMAIL@YOUR.DOMAIN.com                                    
#SBATCH --mail-type=END                           
#SBATCH --mail-type=FAIL                          
#SBATCH -o "stdout.%x.%j.%N"                      
#SBATCH -e "stderr.%x.%j.%N"   
```

Notice the shbang line at the top denoting that this is a bash script and the following commented information below the shbang line. This information is set to a single standard node on Atlas with 48 cores and 7 days of computation time. If you want more information on available nodes an easy way to generate jobs for atlas [look here](https://www.hpc.msstate.edu/computing/atlas/#ondemand). There may be a point where this script requires greater computational resources and you can look a the previous link to assess your needs. Make sure that when you implement this code that you set the email to your own email!

## Meat and Potatoes of the Script

Here are all the options which require modification in the script:

```bash
# Define necessary variables
work_dir="/project/90daydata/guedira_seq_map/atlas_gbs_pipeline_test"
study_name="test_run"
keyfile="/project/guedira_seq_map/zjwinn_working_directory/HPC-GBS-Pipeline/GAM2023_yr-pop_keyfile.txt"
fastq_dir="/project/90daydata/guedira_seq_map/Wheat"
ref_file="/project/90daydata/guedira_seq_map/Refseq2.1_IWGSC/assembly/iwgsc_refseqv2.1_assembly.fa"
enzymes="PstI-MspI"
taglength="85"
ram="300g"
taxamiss="0.85"
taxahet="0.3"
maxdep="100"
snpmiss="0.50"
snphet="0.1"
maf="0.05"
removechr="UNKNOWN"
ncores="40"
disc_and_prod="true"
filter="true"
impute="true"

# Report current working directory
script_dir="/home/zachary.winn/guedira_seq_map/zjwinn_working_directory/HPC-GBS-Pipeline"
```
### work_dir

This option represents the directory where you want all relevant outputs to be written to. Make sure your allotment in those directories are sufficent because discoveries can be quite large.

### study_name

This option denotes what you want the study name to be. ***BE SURE YOU AVOID SPACES!***

### keyfile

This is a tab delimited file with relevant information to the TASSEL pipeline. If you are unfamiliar with this, please refere [here](https://www.maizegenetics.net/copy-of-tassel).

### fastq_dir

This is the directory where all neccesary FASTQ files are located.

### ref_file

This is the true path to a reference genome (e.g., [Chinese Spring Wheat RefSeq Version 2.1](https://urgi.versailles.inrae.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v2.1/)) and it must be indexed with [bwa](https://bio-bwa.sourceforge.net/). 

### enzymes="PstI-MspI"

### taglength="85"

### ram="300g"

### taxamiss="0.85"

### taxahet="0.3"

### maxdep="100"

### snpmiss="0.50"

### snphet="0.1"

### maf="0.05"

### removechr="UNKNOWN"

### ncores="40"

### disc_and_prod="true"

### filter="true"

### impute="true"

### script_dir
