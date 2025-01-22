#!/bin/bash
#SBATCH --job-name="tassel_pl"                   
#SBATCH --qos=normal
#SBATCH -p atlas                                 
#SBATCH -A guedira_seq_map                       
#SBATCH -N 1                                      
#SBATCH -n 48                                     
#SBATCH -t 7-00:00:00                             
#SBATCH --mail-user=zjwinn@ncsu.edu                                      
#SBATCH --mail-type=END                           
#SBATCH --mail-type=FAIL                          
#SBATCH -o "stdout.%x.%j.%N"                      
#SBATCH -e "stderr.%x.%j.%N"                      

# Load necessary modules
module load beagle 
module load bwa
module load samtools
module load bcftools
module load vcftools
module load sqlite
module load bowtie2
module load bzip2

# Directory and file names
directory="tassel-5-standalone"
run_script="tassel-5-standalone/run_pipeline.pl"
repo_url="https://bitbucket.org/tasseladmin/tassel-5-standalone.git"

# Check if the directory exists in the current working directory
if [ ! -d "$directory" ]; then
    echo "$directory not found. Cloning repository..."
    git clone "$repo_url"
else
    echo "$directory already exists."
fi

# Define tassel in the user's path
export TASSEL_PL=$(realpath "$run_script")

# Define necessary variables
work_dir="/project/90daydata/guedira_seq_map/atlas_gbs_pipeline_test"
study_name="test_run"
keyfile="/project/guedira_seq_map/zjwinn_working_directory/HPC-GBS-Pipeline/softWheat_test_discovery_22991.txt"
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

# Echo
echo "The current directory of this bash script is $script_dir"

# Check for discovery and production
if [ $disc_and_prod = "true" ]; then

    # Source the bash script
    source ./tassel_disc_plus_prod_bwa.sh --workdir="$work_dir" \
                                          --study="$study_name" \
                                     	  --keyfile="$keyfile" \
                                     	  --fastq="$fastq_dir" \
                                     	  --ref="$ref_file" \
                                     	  --enzymes="$enzymes" \
                                      	  --taglength="$taglength" \
                                      	  --ram="$ram" \
                                     	  --run-prod

fi

# Change directory back to this directory
cd $script_dir

# Check for filter
if [ $filter = "true" ]; then
    
    # Make directory
    mkdir $work_dir/filt_VCF

    # Run filtering
    source ./vcf_filter_parallel.sh --workdir="${WKDIR}" \
                                    --vcfin=$work_dir/raw_VCF/${study_name}_production.vcf.gz \
                                    --vcfout=$work_dir/filt_VCF/${study_name}_production_filt.vcf.gz \
                                    --taxamiss=$taxamiss \
                                    --taxahet=$taxahet \
                                    --maxdep=$maxdep \
                                    --snpmiss=$snpmiss \
                                    --snphet=$snphet \
                                    --maf=$maf \
                                    --removechr=$removechr \
                                    --ncores=$ncores

fi

# Change directory back to this directory
cd $script_dir

# Check for impute
if [ $impute = "true" ]; then

    # Make directory
    mkdir $work_dir/imp_VCF

    # Run beagle
    beagle gt=$work_dir/filt_VCF/${study_name}_production_filt.vcf.gz \
           out=$work_dir/imp_VCF/${study_name}_production_filt_imp \
           map=./SynOp_RIL906_v1.0_GBS_monotonic.map \
           nthreads=$ncores \
           window=205 

fi
