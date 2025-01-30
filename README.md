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
work_dir="/project/90daydata/guedira_seq_map/atlas_gbs_pipeline_test" # Where the pipeline directory will be made
study_name="test_run" # The name of the project which will be named 
database="/home/zachary.winn/guedira_seq_map_90daydata/atlas_gbs_pipeline_test/database/test_run.db" # This is a path to an existing database, leave blank if discovering and producing in the same run
keyfile="/project/guedira_seq_map/zjwinn_working_directory/HPC-GBS-Pipeline/softWheat_test_discovery_22991.txt" # This is the path to the keyfile
fastq_dir="/project/90daydata/guedira_seq_map/Wheat" # This is the path to the fastq files
ref_file="/project/90daydata/guedira_seq_map/Refseq2.1_IWGSC/assembly/iwgsc_refseqv2.1_assembly.fa" # This is the path to the BWA indexed genome
enzymes="PstI-MspI" # This is enzymes used for the pipeline (usually leave alone)
taglength="85" # This is the minimum length of a read (usually leave alone)
ram="300g" # This is the amount of RAM (edit if you need more and change nodes)
taxamiss="0.85" # Maximum missing data in a line
taxahet="0.3" # Maximum heterozygosity in a line
maxdep="100" # Maximum depth per SNP
snpmiss="0.50" # Maximum missing data per SNP
snphet="0.1" # Maximum heterozygosity per SNP
maf="0.05" # Minimum allele frequency per SNP
removechr="UNKNOWN" # Chromosome to remove (usually leave alone)
ncores="40" # Number of cores for parallele processing (edit if you need more and change nodes)
discovery="false" # Logical to perform discovery
production="false" # Logical to perform production
filter="false" # Logical to perform filter
impute="true" # Logical to do imputation

# Report current working directory
# Note: make sure to set this every time to the scripts working directory!
script_dir="/home/zachary.winn/guedira_seq_map/zjwinn_working_directory/HPC-GBS-Pipeline"
```
### work_dir

This option represents the directory where you want all relevant outputs to be written to. Make sure your allotment in those directories are sufficent because discoveries can be quite large.

### study_name

This option denotes what you want the study name to be. ***BE SURE YOU AVOID SPACES!***

### database

This is the realpath to an exisiting database (.db) of fastq SNPs. Leave this blank if you are running the whole pipeline.

### keyfile

This is a tab delimited file with relevant information to the TASSEL pipeline. If you are unfamiliar with this, please refere [here](https://www.maizegenetics.net/copy-of-tassel).

### fastq_dir

This is the directory where all neccesary FASTQ files are located.

### ref_file

This is the true path to a reference genome (e.g., [Chinese Spring Wheat RefSeq Version 2.1](https://urgi.versailles.inrae.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v2.1/)) and it must be indexed with [bwa](https://bio-bwa.sourceforge.net/). 

### enzymes="PstI-MspI"

This is the type of enzyme used to digest the DNA before ligation to adaptors in sequencing. This enzyme is specific to GBS for wheat. Please make sure this restriction enzyme is correct.

### taglength="85"

This is the minimum length of a read to be considered for alignment. The default for TASSEL is 65, however the regional genotyping lab has elevated this to improve read alignment.

### ram="300g"

This is an added option to specify the maximum RAM available on the computing node. Make sure to alter to fit node specifications and potential RAM requirments.

### taxamiss="0.85"

The total missing data allowed for a single taxa (sample). 

### taxahet="0.3"

The total heterozygosity allowed for a single taxa (sample).

### maxdep="100"

The maximum read depth for a SNP in the filtered VCF.

### snpmiss="0.50"

The maximum amount of missing data allowed for a SNP in a filtered VCF.

### snphet="0.1"

The maximum amount of heterozygosity allowed for a SNP in a filtered VCF.

### maf="0.05"

The minimum allele frequency per SNP in a filtered VCF.

### removechr="UNKNOWN"

The chromosome to filter out of a filtered VCF.

### ncores="40"

The total number of cores available on the computing node.

### discovery="true"

A logical argument to perform discovery. Set to true to perform discovery and production.

## production="true"

A logical argument to perform production. Set to true to perform discovery and production.

### filter="true"

A logical argument to perform filtration. Set to true to perform filtration on a produced VCF.

### impute="true"

A logical argument to perform imputation. Set to true to perform imputation of a filtered VCF.

### script_dir

This is the directory of the git repository. Set this manually every time you run this pipeline!
