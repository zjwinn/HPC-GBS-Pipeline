#!/bin/bash
set -e

#### TASSEL 5 GBS production only
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
## 
## This script runs the TASSEL-GBS SNP-calling pipeline using an input set of
## Fastq files, a keyfile identifying samples from the pool of Fastq files, and
## a database created by TASSEL during the SNP discovery phase.
##
## This version of the script is designed to take command line arguments as input.
## Run the script without any arguments (or with -h or --help) to see the help
## message.
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
opts=$(getopt -o w:d:s:k:f:r:x::e::t::bh --long workdir:,database:,study:,keyfile:,fastq:,ref:,ram::,enzymes::,taglength::,bcf,help -n 'option-parser' -- "$@")
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
--database=<file>   Path to SNP discovery database, either relative to workdir
                    or absolute
--study=<string>    Name of the study (to be added to output files)
--keyfile=<file>    Path to keyfile, either relative to workdir or absolute
--fastq=<dir>       Directory for FASTQ files, either relative to workdir or
                    absolute.

OPTIONAL ARGUMENTS:

--ref=<file>            Path to the uncompressed reference genome FASTA file. 
                        used in the SNP discovery pipeline. This is used to
                        create contig lines in the output VCF (highly recommended)
--enzymes='PstI-MspI'   Restriction enzymes used to generated GBS library
--taglength=64          The maximum length of tags used in the Tassel database.
                        Note that this should be no more than the length of the
                        sequencing reads minus the length of the longest barcode.
                        In practice it should be slightly shorter than this to
                        trim off low-quality bases from the end of Illumina reads.
                        The default 64bp is computationally efficient but may
                        not be ideal for low-diversity panels/populations/organisms.
--bcf                   Output a BCF file, instead of the default bgzipped VCF file


HELP ARGUMENTS:
--help or -h    Displays this message."

echo -e \\n"USAGE EXAMPLE: 
${bold}./$script --workdir=~/my_gbs_run \\
                    --database=~/path/to/discovery/db \\
                    --study=GBS_run01 \\
                    --keyfile=GBS_run01_keyfile.txt \\
                    --fastq=~/path/to/fastq/files/directory \\
                    --ref=~/path/to/ref_genome.fasta \\
                    --enzymes=PstI-MspI \\
                    --taglength=64 \\
                    --ram=300g ${norm}"

echo "
DEPENDENCIES:
    1)  TASSEL 5, with the absolute path to run_pipeline.pl defined as a system
        variable named TASSEL_PL
    2)  samtools - used to create contig lines in output VCF
    3)  bcftools - used to create and index final .vcf.gz or .bcf output file
    4)  VCFtools (optional) - used to calculate depth statistics at end. If not
        present, the calculation of depth statistics will be skipped.

OUTPUT: 
    The script creates the directory raw_VCF/ within the specified working
    directory, and places all output there.

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
TAG_LENGTH=64
RG="NA"
ext="vcf"

## Parse out command-line arguments
while true; do
    case "$1" in
        -w | --workdir )
            WD="$2" 
            shift 2
            ;;
        -d | --database ) 
            DB="$2" 
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
        -x | --ram )
            RAM="$2"
            shift 2
            ;;
        -b | --bcf ) 
            ext="bcf" 
            shift
            ;;
        -h | --help ) 
            help
            ;;
        -- ) 
            shift 
            break
            ;;
        * ) echo "Internal error!" break;;
    esac
done


#### Run TASSEL-GBS ###########################################################

## Change into the working directory
mkdir -p $WD
cd $WD

## Convert paths to absolute
KF="$(realpath ${KF})"
DB="$(realpath ${DB})"
FASTQ="$(realpath ${FASTQ})"
if [[ $RG != "NA" ]]; then
    RG="$(realpath ${RG})"
fi

#### Check for existence of user-supplied files/dirs ####
if [[ ! -f "$KF" ]]; then
    echo "The specified keyfile does not exist!"
    exit 1;
fi

if [[ ! -f "$DB" ]]; then
    echo "The specified discovery database does not exist!"
    exit 1;
fi

if [[ ! -d "$FASTQ" ]]; then
    echo "The specified fastq directory does not exist!"
    exit 1;
fi

if [[ ! -f "$RG" ]]; then
    echo "WARNING - The specified reference genome fasta file does not exist!"
fi

## One annoying idiosyncracy - the FASTQ directory must have a trailing slash
if [[ "${FASTQ: -1}" != "/" ]]; then
    FASTQ=${FASTQ}/
fi


## Record System and Run Parameters
echo "OS version:" > system_run_info.txt
cat /etc/os-release >> system_run_info.txt
echo -e "\nReference genome & md5 sum:" >> system_run_info.txt
md5sum "${RG}" >> system_run_info.txt
echo -e "\nTASSEL version:" >> system_run_info.txt
$TASSEL_PL -h 2>&1 | grep "Tassel Version" >> system_run_info.txt
echo -e "\nBCFTools/HTSlib version:" >> system_run_info.txt
bcftools --version >> system_run_info.txt
echo -e "\nVCFTools version:" >> system_run_info.txt
vcftools --version >> system_run_info.txt
if [[ -f "${beajar}" ]]; then
    echo -e "\nBeagle version:" >> system_run_info.txt
    echo $(basename "${beajar}") >> system_run_info.txt
fi


## Create file of contig lines using samtools and awk
if [[ -f "$RG" ]]; then
    if [[ ! -f "$RG".fai ]]; then
        samtools faidx $RG
    fi
    awk 'BEGIN {OFS = ""}{print "##contig=<ID=", $1, ",length=", $2, ">"}' "$RG".fai > contig_lines.txt
fi

## Create TASSEL output directory
mkdir -p ./logs ./raw_VCF

## Create the name of the output files
OUTFILE="raw_VCF/${STUDY}_production.vcf"
if [[ "$ext" == "vcf" ]]; then
    FINAL_FILE="${OUTFILE}.gz"
else
    FINAL_FILE=$(echo "$OUTFILE" | sed 's/vcf$/bcf/')
fi

## Run Production
echo "Production" > ./logs/production.log
date >> ./logs/production.log
$TASSEL_PL -Xms25g -Xmx$RAM -fork1 -ProductionSNPCallerPluginV2 -db $DB -e $E -i $FASTQ -k $KF -kmerLength $TAG_LENGTH -o $OUTFILE -endPlugin -runfork1 > ./logs/ProductionSNPCallerPluginV2.log


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

## Add contig info if present
## Unfortunately Tassel seems to convert any lowercase letters in chrom/contig
## IDs to uppercase, which has caused unexpected headaches
if [[ -f contig_lines.txt ]]; then
    cat contig_lines.txt >> raw_VCF/header.vcf
    sed -e 's/##contig=<ID=.*,/\U&/' raw_VCF/header.vcf
    sed -e 's/##CONTIG/##contig/' raw_VCF/header.vcf
fi

## Update VCF with new header lines
grep -v "^##" $OUTFILE >> raw_VCF/header.vcf
mv raw_VCF/header.vcf $OUTFILE


#### Compress/Index VCF #######################################################

## Compress VCF file to either .vcf.gz or .bcf file
if [[ "$ext" == "vcf" ]]; then
    bcftools view "$OUTFILE" -Oz -o "$FINAL_FILE"
else
    bcftools view "$OUTFILE" -Ob -o "$FINAL_FILE"
fi
bcftools index -f "$FINAL_FILE"


#### Calculate summary stats ##################################################

echo "Genotype summary" >> ./logs/production.log
date >> ./logs/discovery.log
$TASSEL_PL -Xmx$RAM \
        -vcf "$OUTFILE" \
        -GenotypeSummaryPlugin \
        -endPlugin \
        -export raw_VCF/summary


#### Calculate depth statistics ################################################

vcftools --gzvcf "$OUTFILE" --site-mean-depth --out raw_VCF/raw_VCF

rm "$OUTFILE"
echo "Script finished at:" >> ./logs/production.log
date >> ./logs/production.log

# Commenting this out so that you can pipeline this script
#exit 0;
