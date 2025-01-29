#!/bin/bash
set -e
LC_ALL=C

## Filter VCF files
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
## 
## Rakin Rouf
## rakinr57@gmail.com
##
## This script is my first attempt at implementing command-line control of a
## bash script. See the help documentation (run the script without any arguments,
## or with -h or --help for details on input/output)
##
## Dependencies:
##   1) VCFTools in user's PATH as "vcftools"
##   2) bcftools in user's PATH as "bcftools".
##   2) PLINK 1.9 (or PLINK 2 when it becomes beta) in user's PATH as "plink2"
##        this is only necessary if filtering by LD
##   3) TASSEL - this is required only to generate summary statistics of
##        filtered data at end. Path to run_pipeline.pl must be modified below 
##        the script will still output filtered VCF if TASSEL fails to
##        start. Tassel run_pipeline.pl should be pointed to by an exported
##        alias, $TASSEL_PL
##
##
## Notes from Zachary Winn
## zwinn@outlook.com
## https://github.com/zjwinn
##
## I am not changing much in this script, I am commenting out exit 0 at the end
###############################################################################

## Path to start_tassel.pl
#tassel="/TASSEL/tassel5v44/run_pipeline.pl"

## Set Script Name variable
script="$(basename ${BASH_SOURCE[0]})"

## Template getopt string below
## Short options specified with -o. option name followed by nothing indicates
## no argument (i.e. flag) one colon indicates required argument, two colons
## indicate optional argument
##
## Long arguments specified with --long. Options are comma separated. Same
## options syntax as for short options
opts=$(getopt -o w:i:o:s::m::n::a::t::e::b::r::l::d::p::c::h --long workdir:,vcfin:,vcfout:,taxasub::,maf::,snpmiss::,snphet::,taxamiss::,taxahet::,minbp::,removechr::,maxld::,mindep::,maxdep::,ncores::,help -n 'option-parser' -- "$@")
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

--vcfin=<file>      Input VCF file path, relative to workdir or absolute path 
                    (can be gzipped)
--vcfout=<file>     Output VCF file path, relative to workdir or absolute path 
                    (output will always be gzipped, or bgzipped and indexed if
                    bgzip and tabix are installed)

OPTIONAL ARGUMENTS - Directories:

--workdir=<dir>    Working directory

OPTIONAL ARGUMENTS - SNP filtering:

NOTE:   Heterozygosity proportions below calculated with respect to NON-MISSING 
        genotype calls

--maf=0         Minimum minor allele frequency. Default value implements no
                filtering
--snpmiss=1     Maximum SNPwise missing data proportion. Default value 
                implements no filtering.
--snphet=1      Maximum SNPwise heterozygous proportion. Default value 
                implements no filtering.
--minbp=0       Minimum physical distance between any two SNPs. Default value 
                implements no filtering.
--mindep=0      Minimum average SNP sequence depth. Default value 
                implements no filtering.
--maxdep=1e6    Maximum average SNP sequence depth. Default value of 1e6 should
                not implement any filtering
--maxld=NA      Maximum pairwise LD between SNPs on each chromosome. 
                Default value implements no filtering.
--removechr=NA   Remove a specific chromosome - default retains all chromosomes

OPTIONAL ARGUMENTS - Taxa filtering:
--taxasub=NA    Text file indicating taxa to manually keep, listed one 
                per line. Path can be relative to workdir or absolute.
                Default value performs no subsetting.
--taxamiss=1    Maximum Taxawise missing data proportion. Default value 
                implements no filtering.
--taxahet=1     Maximum Taxawise heterozygous proportion. Default value
                implements no filtering.
                
OPTIONAL ARGUMENTS - Number of cores                
--ncores=2      Number of cores to use.

HELP ARGUMENTS:
--help or -h    Displays this message."
    
echo -e \\n"USAGE EXAMPLE: 
${bold}./$script --vcfin=~/geno.vcf --vcfout=~/geno_out.vcf.gz --maf=0.05 ${norm}"

echo "
DEPENDENCIES:
    1) VCFTools in user's PATH as 'vcftools'
    2) bcftools in user's PATH as 'bcftools'
    2) PLINK 1.9 (or PLINK 2 when it becomes beta) in user's PATH as 'plink2'
        only required for filtering by LD
    3) TASSEL for generating filtered file statistics - path to run_pipeline.pl
        is set within script. Script will still create filtered VCF without
        using TASSEL
"
exit 1;
}

## If no arguments supplied, print help message
if [[ $1 == "--" ]]; then help; fi

## Set initial values for optional arguments
vcfin="init"
vcfout="init"
workdir="NA"
taxsub="NA"
maf=0
snpmiss=1
snphet=1
taxmiss=1
taxhet=1
mindep=0
maxdep=1e6
minbp=0
removechr="NA"
maxld="NA"
N=2

## Parse out command-line arguments
while true; do
    case "$1" in
    -w | --workdir ) workdir="$2"; shift 2;;
    -i | --vcfin ) vcfin="$2"; shift 2;;
	-o | --vcfout ) vcfout="$2"; shift 2;;
    -s | --taxasub ) taxsub="$2"; shift 2;;
	-m | --maf ) maf="$2"; shift 2;;
	-n | --snpmiss ) snpmiss="$2"; shift 2;;
    -a | --snphet ) snphet="$2"; shift 2;;
	-t | --taxamiss ) taxmiss="$2"; shift 2;;
    -e | --taxahet ) taxhet="$2"; shift 2;;
    -b | --minbp ) minbp="$2"; shift 2;;
    -r | --removechr ) removechr="$2"; shift 2;;
    -l | --maxld ) maxld="$2"; shift 2;;
    -d | --mindep ) mindep="$2"; shift 2;;
    -p | --maxdep ) maxdep="$2"; shift 2;;
    -c | --ncores ) N="$2"; shift 2;;
    -h | --help ) help;;
    -- ) shift; break;;
    * ) echo "Internal error!" break;;
    esac
done

## If user-requested cores exceed the maximum, set to the number of cores to (maximum cores - 1)
max_cores=$(nproc)
if [ $N -gt $max_cores ]; then
	let "N = $max_cores - 1"
	echo "maximum cores exceeded, setting number of cores to $N"
fi

## Check input and output files were provided
if [[ "$vcfin" == "init" ]] || [[ "$vcfout" == "init" ]]; then
    echo "--vcfin and --vcfout both require arguments. Type './vcf_filter_command_opts.sh -h' for help"
    exit 1;
fi

## VCFTools filters SNPs by min non-missing data, which I find unintuitive
## compared to max missing data - so we need to find the opposite ratio
snpmissinv="$(echo "1 - $snpmiss" | bc)"

## Change to working directory if one was supplied
if [[ "$workdir" != "NA" ]]; then cd "$workdir"; fi

## Remove .gz extension from output filename if necessary
outext="${vcfout##*.}"
if [[ "$outext" == "gz" ]]; then
    vcfout="${vcfout::-3}"
fi 

## Create output directory if necessary and temp directory within it
outdir="$(dirname ${vcfout})"
mkdir -p "$outdir"

## Convert paths to absolute
outdir="$(realpath ${outdir})"
vcfin="$(realpath ${vcfin})"
vcfout="$(realpath ${vcfout})"
if [[ "$taxsub" != "NA" ]]; then taxsub="$(realpath ${taxsub})"; fi

## Create temp dir in output dir
tempdir="$(mktemp -d -p "$outdir")"
if [[ ! "$tempdir" || ! -d "$tempdir" ]]; then
    echo "Could not create temporary directory"
    exit 1;
fi

## Write filtering parameters to table

echo -e "parameter\tvalue
vcfin\t$vcfin
vcfout\t$vcfout
taxa_sub_list\t$taxsub
taxa_max_miss\t$taxmiss
taxa_max_het\t$taxhet
min_maf\t$maf
snp_max_miss\t$snpmiss
snp_max_het\t$snphet
snp_min_dep\t$mindep
snp_max_dep\t$maxdep
min_bp_space\t$minbp
remove_chromosome\t$removechr
max_ld\t$maxld" > "${outdir}/filtering_params.tsv"

cd "$tempdir"
#echo "Changed to working directory: ${PWD}"

## Extract extension from (possibly gzipped) input VCF
ext="${vcfin##*.}"

#### FILTERING ROUND 1 - Filter taxa, clean up SNPs ###########################

## Generate lists of Taxa to keep by missingness and heterozygosity
echo "Calculating taxa-wise missingness and heterozygosity"
if [ "$ext" == "gz" ]; then 
	vcftools --gzvcf "$vcfin" --missing-indv
	vcftools --gzvcf "$vcfin" --het
else
	vcftools --vcf "$vcfin" --missing-indv
	vcftools --vcf "$vcfin" --het
fi

## Find taxa that pass the missing data and heterozygosity filters
awk -F"\t" -v taxmiss="$taxmiss" '$5 <= taxmiss { print $1 }' out.imiss \
	| sort > tax_miss_keep.txt
tail -n +2 out.het \
	| awk -F"\t" -v taxhet="$taxhet" '($4-$2)/$4 <= taxhet { print $1 }' \
	| sort > tax_het_keep.txt

## Intersect the two lists above
comm -12 tax_miss_keep.txt tax_het_keep.txt > taxa_keep.txt

## If taxsub was provided, intersect it with taxa_keep.txt
if [[ "$taxsub" != "NA" ]]; then
    comm -12 <(sort $taxsub) <(sort taxa_keep.txt) > taxa_temp.txt
    mv taxa_temp.txt taxa_keep.txt
fi

## Perform the initial filtering
echo "Filtering by taxa"
if [[ "$ext" == "gz" ]]; then      
    vcftools --gzvcf "$vcfin" \
                --keep taxa_keep.txt \
                --min-alleles 2 \
                --max-alleles 2 \
                --remove-indels \
                --recode \
                --out round1
else
    vcftools --vcf "$vcfin" \
                --keep taxa_keep.txt \
                --min-alleles 2 \
                --max-alleles 2 \
                --remove-indels \
                --recode \
                --out round1
fi


#### FILTERING ROUND 2 - Filter SNPs ##########################################

## Create temp directory to hold chromosome files
mkdir ./temp_storage

## Function for getting the number of het calls and filtering
GetHetCalls() {

	chrom=$1
	
	if	[[ $snphet == 1 ]]; then
		grep "^$chrom" round1.recode.vcf | cut -f 3 > temp_storage/snps_keep_$chrom.txt
	else
		echo "Calculating number of non-missing calls for $chrom"
		## Calculate number of non-missing calls per SNP
		grep -E "^#|^$chrom" round1.recode.vcf | bcftools query -f '[\S%CHROM\_%POS\n]' -i 'GT!~".\."' | sort | uniq -c > temp_storage/non_miss_$chrom.txt

		echo "Calculating number of het calls for $chrom"
		## Calculate number of het calls per SNP
		grep -E "^#|^$chrom" round1.recode.vcf | bcftools query -f '[\S%CHROM\_%POS\n]' -i 'GT="het"' | sort | uniq -c > temp_storage/het_$chrom.txt
	
        ## Join the two files together, inserting 0 for missing vals (should only be
        ## values from the het file that have missing)
        ## Format is: SNP_name het_count non-missing_count
		join -a 2 -e 0 -1 2 -2 2 -o 2.2,1.1,2.1 temp_storage/het_$chrom.txt temp_storage/non_miss_$chrom.txt > temp_storage/joined_$chrom.txt
		
		## Find SNPs passing the het. threshold
		awk -v het="$snphet" '$2/$3 <= het {print $1}' temp_storage/joined_$chrom.txt > temp_storage/snps_keep_$chrom.txt
	fi
	
	## Implement filtering of SNPs
	if [[ "$removechr" == "NA" ]]; then
		vcftools --vcf round1.recode.vcf \
                    --snps temp_storage/snps_keep_$chrom.txt \
                    --maf $maf \
                    --max-missing $snpmissinv \
                    --min-meanDP $mindep \
                    --max-meanDP $maxdep \
                    --thin $minbp \
                    --recode \
                    --out temp_storage/round2_$chrom
	else
		vcftools --vcf round1.recode.vcf \
                    --snps temp_storage/snps_keep_$chrom.txt \
                    --maf $maf \
                    --max-missing $snpmissinv \
                    --min-meanDP $mindep \
                    --max-meanDP $maxdep \
                    --thin $minbp \
                    --not-chr $removechr \
                    --recode \
                    --out temp_storage/round2_$chrom
	fi
}

## Get list of chromosome parts
chrom_parts=$(grep -v "^#" round1.recode.vcf | cut -f 1 | uniq)

echo "Calculating SNP-wise heterozygosity"

## Parallel loop for filtering and getting het calls
for chrom in $chrom_parts
do
	test "$(jobs | wc -l)" -ge $N && wait -n || true
	GetHetCalls "$chrom" &
done
wait


## Get header from a one of the chromosome files
head_file=$(ls temp_storage/round2*.vcf | head -n1)
grep "^#" $head_file > temp_storage/header.txt

## Concat all the seperate chromosome files into round2.recode.vcf
bcftools concat temp_storage/round2*.recode.vcf | bcftools sort -T "$tempdir" | grep -v "^#" | cat temp_storage/header.txt - | sed '1 a ##FILTER=<ID=PASS,Description="All filters passed">' > round2.recode.vcf


#### FILTERING ROUND 3 (Optional) - Filter by LD ##############################
if [[ "$maxld" == "NA" ]]; then
    mv round2.recode.vcf round3.recode.vcf
else
	echo "Filtering by LD"
    plink2 --vcf round2.recode.vcf \
            --id-delim ^ \
            --allow-extra-chr \
            --indep-pairwise 50000 100 $maxld \
            --out thinning

    vcftools --vcf round2.recode.vcf \
            --snps thinning.prune.in \
            --recode \
            --out round3
fi


#### PROCESS OUTPUTS ##########################################################

## Check if bcftools installed, if not use regular gzip
if command -v bcftools; then
	if command -v pbgzip; then
		pbgzip -c -t 0 round3.recode.vcf > "${vcfout}.gz"
	else
        bcftools view -o "${vcfout}.gz" -O z --no-version round3.recode.vcf
    fi
    bcftools index "${vcfout}.gz"
else
    gzip -c round3.recode.vcf > "${vcfout}.gz"
fi

## Remove temp. directory
cd "$outdir"
rm -rf "$tempdir"
echo "Removed temp directory $tempdir"

## Attempt to generate summary statistics using TASSEL
echo "Generating summary statistics with TASSEL"
$TASSEL_PL -vcf "${vcfout}.gz" \
        -GenotypeSummaryPlugin \
        -endPlugin \
        -export summary

#### Calculate depth statistics ################################################
vcftools --gzvcf ${vcfout}.gz --site-mean-depth --out ./filtered_VCF        

# Comment out the last exit so that this can be pipelined together with other scripts
#exit 0;
