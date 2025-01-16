#!/bin/bash
set -e

#### TASSEL 5 GBS discovery plus production BWA - command line input
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
## 
## This script runs the TASSEL-GBS SNP-calling pipeline using an input set of
## Fastq files, a reference genome, and a keyfile identifying samples from the
## pool of Fastq files.
##
## This script always performs SNP discovery - it optionally also runs production
## on the same set of lines supplied in the keyfile if the --run-prod flag is supplied.
## Running a separate large discovery phase and subsequent production runs is
## often advantageous when working with large diversity panels. For mapping
## populations, we often run discovery and production, using just the lines of the
## population.
##
## This version of the script is designed to take command line arguments as input.
## Run the script without any arguments (or with -h or --help) to see the help
## message.
##
## Zachary Winn
## zwinn@outlook.com
## https://github.com/zjwinn
##
## Edit notes from 1/15/2024
## I have added a porition that requires the user to define the amount of RAM
## So that this pipeline can be used on Atlas or Ceres from SCINet. Also,
## I have written a partner SLURM script to submit this pipeline as a SLURM 
## job. This should allow the program to be used on the HPC in a pinch.
###############################################################################

## Set Script Name variable
script="$(basename ${BASH_SOURCE[0]})"

## Template getopt string below
## Short options specified with -o. option name followed by nothing indicates
## no argument (i.e. flag) one colon indicates required argument, two colons
## indicate optional argument
##
## Long arguments specified with --long. Options are comma separated. Same
## options syntax as for short options
opts=$(getopt -o w:s:k:f:r:e:x::ipt::m::h --long workdir:,study:,keyfile:,fastq:,ref:,enzymes:,ram::,inclref,run-prod,taglength::,minq::,help -n 'option-parser' -- "$@")
eval set -- "$opts"

## Set fonts used for Help.
norm="$(tput sgr0)"
bold="$(tput bold)"
rev="$(tput smso)"

## Help function
function help {
echo -e \\n"Help documentation for ${bold}${script}.${norm}"\\n
echo "
REQUIRED ARGUMENTS:

--workdir=<dir>     Working directory
--study=<string>    Name of the study (to be added to output files)
--keyfile=<file>    Path to keyfile, either relative to workdir or absolute
--fastq=<dir>       Directory for FASTQ files, either relative to workdir or
                    absolute.
--ref=<file>        Path to an uncompressed reference genome FASTA file. Note that
                    this file must be located in the same directory as its
                    corresponding BWA index file. 
--ram=<int>'g'      This is an integer that represents the amount of RAM available in gigabytes 
                    followed by the string 'g'.

OPTIONAL ARGUMENTS:

--run-prod              Flag indicating whether to run production using the
                        supplied keyfile. Otherwise only SNP discovery is run.
--enzymes='PstI-MspI'   Restriction enzymes used to generated GBS library
--taglength=64          The maximum length of tags used in the Tassel database.
                        Note that this should be no more than the length of the
                        sequencing reads minus the length of the longest barcode.
                        In practice it should be slightly shorter than this to
                        trim off low-quality bases from the end of Illumina reads.
                        The default 64bp is computationally efficient but may
                        not be ideal for low-diversity panels/populations/organisms.
--minq=0                The minimum quality score within a tag. Trimmed tags where
                        any base (including barcodes) falls below this threshold are
                        removed. Default disables quality-based filtering.
--inclref               Flag indicating whether to include a tag
                        from the reference genome when calling SNPs in production

HELP ARGUMENTS:
--help or -h    Displays this message."
    
echo -e \\n"USAGE EXAMPLE: 
${bold}./$script --workdir=~/my_gbs_run \\
                                   --study=GBS_run01 \\
                                   --keyfile=GBS_run01_keyfile.txt \\
                                   --fastq=~/path/to/fastq/files/directory \\
                                   --ref=~/path/to/ref_genome.fasta \\
                                   --enzymes=PstI-MspI \\
                                   --taglength=64 \\
                                   --minq=0 \\
                                   --ram=300g \\
                                   --run-prod \\
                                   --inclref ${norm}"

echo "
DEPENDENCIES:
    1) TASSEL 5, with the absolute path to run_pipeline.pl defined as a system
       variable named TASSEL_PL
    2) samtools installed in user's PATH
    3) bgzip (optional) - used to compress output VCF. If not installed,
       regular gzip is used
    4) bcftools (optional) - used to index output VCF
    5) VCFtools (optional) - used to calculate depth statistics at end. If not
       present, the calculation of depth statistics will be skipped.

OUTPUT: 
    The script creates 4 sub-directories in the specified working directory:
        1) logs/
        2) database/
        3) alignment/
        4) raw_VCF/ (only if running production)

NOTES:
    This script will correct a few problems with the header lines of VCF files
    produced by TASSEL, including adding in contig lines.
"
exit 1;
}

## If no arguments supplied, print help message
if [[ $1 == "--" ]]; then help; fi

## Set initial values for optional arguments
E="PstI-MspI"
INCL_REF="false"
RUN_PROD="false"
TAG_LENGTH=64
MIN_Q=0

## Parse out command-line arguments
while true; do
    case "$1" in
        -w | --workdir )
            WD="$2"
            shift 2
            ;;
        -s | --study ) 
            STUDY="$2"
            shift 2
            ;;
	    -k | --keyfile ) 
            KF="$2"
            shift 2
            ;;
        -f | --fastq ) 
            FASTQ="$2" 
            shift 2
            ;;
	    -r | --ref ) 
            RG="$2" 
            shift 2
            ;;
	    -e | --enzymes )
            E="$2"
            shift 2
            ;;
	    -t | --taglength )
            TAG_LENGTH="$2"
            shift 2
            ;;
        -m | --minq )
            MIN_Q="$2"
            shift 2
            ;;
        -x | --ram )
            RAM="$2"
            shift 2
            ;;
        -i | --inclref ) 
            INCL_REF="true"
            shift
            ;;
        -p | --run-prod )
            RUN_PROD="true"
            shift
            ;;
        -h | --help ) 
            help
            ;;
        -- ) 
            shift
            break
            ;;
        * )
            echo "Internal error!" 
            break
            ;;
    esac
done


#### Run TASSEL-GBS ###########################################################

## Change into the working directory
mkdir -p "${WD}"
cd "${WD}"

## Convert paths to absolute
KF="$(realpath ${KF})"
FASTQ="$(realpath ${FASTQ})"
RG="$(realpath ${RG})"

#### Check for existence of user-supplied files/dirs ####
if [ ! -f "${KF}" ]; then
    echo "The specified keyfile does not exist!"
    exit 1;
fi

if [ ! -d "${FASTQ}" ]; then
    echo "The specified fastq directory does not exist!"
    exit 1;
fi

if [ ! -f "${RG}" ]; then
    echo "The specified reference genome fasta file does not exist!"
    exit 1;
fi

## One annoying idiosyncracy - the FASTQ directory must have a trailing slash
if [[ "${FASTQ: -1}" != "/" ]];then
    FASTQ="${FASTQ}"/
fi


## Record System and Run Parameters
echo "OS version:" > system_run_info.txt
cat /etc/os-release >> system_run_info.txt
echo -e "\nReference genome & md5 sum:" >> system_run_info.txt
md5sum "${RG}" >> system_run_info.txt
echo -e "\nTASSEL version:" >> system_run_info.txt
$TASSEL_PL -h 2>&1 | grep "Tassel Version" >> system_run_info.txt
echo -e "\nBWA version:" >> system_run_info.txt
bwa 2>&1 | grep "Version" >> system_run_info.txt
echo -e "\nBCFTools/HTSlib version:" >> system_run_info.txt
bcftools --version >> system_run_info.txt
echo -e "\nVCFTools version:" >> system_run_info.txt
vcftools --version >> system_run_info.txt
if [[ -f "${beajar}" ]]; then
    echo -e "\nBeagle version:" >> system_run_info.txt
    echo $(basename "${beajar}") >> system_run_info.txt
fi


## Create file of contig lines using samtools and awk
if [ ! -f "${RG}".fai ]; then
    samtools faidx "${RG}"
fi
awk 'BEGIN {OFS = ""}{print "##contig=<ID=", $1, ",length=", $2, ">"}' "${RG}".fai > contig_lines.txt

## Create TASSEL output directories
mkdir -p ./logs ./database ./alignment
if [[ "${RUN_PROD}" == "true" ]]; then
    mkdir ./raw_VCF
fi

## Create the name of the output file
OUTFILE="raw_VCF/${STUDY}_production.vcf"

echo "Identify tags and output to database file" > ./logs/discovery.log
date >> ./logs/discovery.log
$TASSEL_PL -Xms25g -Xmx$RAM -fork1 -GBSSeqToTagDBPlugin -e $E -i $FASTQ -db ./database/${STUDY}.db -k $KF -kmerLength $TAG_LENGTH -mnQS $MIN_Q -c 5 -mxKmerNum 50000000 -deleteOldData true -endPlugin -runfork1 > ./logs/GBSSeqToTagDBPlugin.log

echo "Retrieve distinct tags from database and reformat to .fq" >> ./logs/discovery.log
date >> ./logs/discovery.log
$TASSEL_PL -Xms25g -Xmx$RAM -fork1 -TagExportToFastqPlugin -db database/${STUDY}.db -o alignment/${STUDY}_MasterGBStags.fa.gz -endPlugin -runfork1 > ./logs/TagExportToFastqPlugin.log

echo "Align the tags to the reference genome using BWA" >> ./logs/discovery.log
date >> ./logs/discovery.log
ncores=$(( $(nproc) - 2 ))
bwa mem -t $ncores $RG alignment/${STUDY}_MasterGBStags.fa.gz > alignment/${STUDY}_AlignedMasterTags.sam

echo "Read SAM file to determine potential positions of Tags against Ref Genome" >> ./logs/discovery.log
date >> ./logs/discovery.log
$TASSEL_PL -Xms25g -Xmx$RAM -fork1 -SAMToGBSdbPlugin -i alignment/${STUDY}_AlignedMasterTags.sam -db database/${STUDY}.db -aLen 0 -aProp 0.0 -endPlugin -runfork1 > ./logs/SAMToGBSdbPlugin.log

echo "Call SNPs" >> ./logs/discovery.log
date >> ./logs/discovery.log

## Check for included reference option
if [[ "$INCL_REF" == "true" ]]; then
	$TASSEL_PL -Xms25g -Xmx$RAM -fork1 -DiscoverySNPCallerPluginV2 -db database/${STUDY}.db -mnMAF 0.01 -mnLCov 0.1 -deleteOldData true -ref $RG -endPlugin -runfork1 > ./logs/DiscoverySNPCallerPluginV2.log
else
	$TASSEL_PL -Xms25g -Xmx$RAM  -fork1 -DiscoverySNPCallerPluginV2 -db database/${STUDY}.db -mnMAF 0.01 -mnLCov 0.1 -deleteOldData true -endPlugin -runfork1 > ./logs/DiscoverySNPCallerPluginV2.log
fi

## Check if RUN_PROD is specified as true
if [[ $RUN_PROD == "true" ]]; then
    ## Run production of VCF
    echo "Production" > ./logs/production.log
    date >> ./logs/production.log
    $TASSEL_PL -Xms25g -Xmx$RAM -fork1 -ProductionSNPCallerPluginV2 -db database/${STUDY}.db -e $E -i $FASTQ -k $KF -kmerLength $TAG_LENGTH -o $OUTFILE -endPlugin -runfork1 > ./logs/ProductionSNPCallerPluginV2.log

    #### Format the output VCF ####################################################

    ## There's a bunch of stuff in the TASSEL-generated VCF file's header
    ## lines that can be improved
    grep "^##" $OUTFILE > raw_VCF/header.vcf

    ## If the reference fasta included "chr" before chromosome names, these are
    ## removed by Tassel, so also need to be stripped out of the contig lines 
    ## in the header for bcftools to work
    sed -i 's/ID=chr/ID=/' raw_VCF/header.vcf

    ## Tassel also capitalizes all characters in chrom names. In the case of wheat
    ## if the unaligned contigs are given chromosome designator "Un", this must
    ## be capitalized
    sed -i 's/ID=Un/ID=UN/' raw_VCF/header.vcf

    ## AD and PL FORMAT fields defined incorrectly
    ## Also, TASSEL's description for the PL field contains equal signs and commas, which causes problems with VCFTools
    sed -i 's/AD,Number=./AD,Number=R/g' raw_VCF/header.vcf
    sed -i 's/PL,Number=./PL,Number=G/g' raw_VCF/header.vcf
    sed -i 's/Description="Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic"/Description="Normalized Phred-scaled likelihoods for genotypes as defined in the VCF specification"/g' raw_VCF/header.vcf

    ## AF and NS INFO fields do not appear to actually be present in file - delete
    grep -E -v "AF,Number=|NS,Number=" raw_VCF/header.vcf > raw_VCF/temp.vcf
    mv raw_VCF/temp.vcf raw_VCF/header.vcf

    ## The "QualityScore" INFO field is not defined in header - add
    echo "##INFO=<ID=QualityScore,Number=1,Type=Float,Description=\"Quality Score\">" >> raw_VCF/header.vcf

    ## Add contig info
    ## Unfortunately Tassel seems to convert any lowercase letters in chrom/contig
    ## IDs to uppercase, which has caused unexpected headaches
    cat contig_lines.txt >> raw_VCF/header.vcf
    sed -e 's/##contig=<ID=.*,/\U&/' raw_VCF/header.vcf
    sed -e 's/##CONTIG/##contig/' raw_VCF/header.vcf

    ## Update VCF with new header lines
    grep -v "^##" $OUTFILE >> raw_VCF/header.vcf
    mv raw_VCF/header.vcf $OUTFILE

    #### Compress/Index VCF #######################################################

    ## Compress VCF file
    ## If bcftools is installed, bgzip and index, if not, regular gzip.
    ## bgzip should be installed if bcftools is
    if command -v bcftools; then
	    bgzip -f $OUTFILE
	    bcftools index -f ${OUTFILE}.gz
    else
        gzip -f $OUTFILE
    fi

    #### Calculate summary stats ##################################################

    echo "Genotype summary" >> ./logs/discovery.log
    date >> ./logs/discovery.log
    $TASSEL_PL -Xms25g -Xmx$RAM -vcf ${OUTFILE}.gz -GenotypeSummaryPlugin -endPlugin -export raw_VCF/summary

    #### Calculate depth statistics ################################################

    vcftools --gzvcf ${OUTFILE}.gz --site-mean-depth --out raw_VCF/raw_VCF
fi

echo "Script finished at:" >> ./logs/discovery.log
date >> ./logs/discovery.log


# Comment out the exit so that this step can be pipelined
#exit 0;
