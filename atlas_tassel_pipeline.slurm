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
repo_url="https://bitbucket.org/tasseladmin/tassel-5-standalone.git" # Hopefully this stays consistent... Maybe check if this is not correct!

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

# Echo
echo "The current directory of this bash script is $script_dir"

# Validate required files and directories
if [ ! -f "$keyfile" ]; then
    echo "Error: Keyfile not found at $keyfile. Exiting."
    exit 1
fi
if [ ! -f "$ref_file" ]; then
    echo "Error: Reference file not found at $ref_file. Exiting."
    exit 1
fi
if [ ! -d "$fastq_dir" ]; then
    echo "Error: FASTQ directory not found at $fastq_dir. Exiting."
    exit 1
fi
if [ ! -d "$work_dir" ]; then
    echo "Error: working directory not found at $work_dir. Exiting"
    exit 1
fi

# Check for discovery and production
if [ $discovery = "true" ]; then

    # Run discovery
    source ./tassel_disc_plus_prod_bwa.sh --workdir="$work_dir" \
                                            --study="$study_name" \
                                            --keyfile="$keyfile" \
                                            --fastq="$fastq_dir" \
                                            --ref="$ref_file" \
                                            --enzymes="$enzymes" \
                                            --taglength="$taglength" \
                                            --ram="$ram"

    # Set discovery location post discovery creation
    database="$work_dir/database/$study_name.db"

fi

# Change directory back to this directory
cd $script_dir

if [ $production = "true" ]; then

    # Validate database existence
    if [ ! -f "$database" ]; then
        echo "Error: Database not found at $database. Exiting."
        exit 1
    fi

    # Run production
    source ./tassel_production.sh --workdir="$work_dir" \
                                    --database="$database" \
                                    --study="$study_name" \
                                    --keyfile="$keyfile" \
                                    --fastq="$fastq_dir" \
                                    --ref="$ref_file" \
                                    --enzymes="$enzymes" \
                                    --taglength="$taglength" \
                                    --ram="$ram"

fi

# Change directory back to this directory
cd $script_dir

# Check for filter
if [ $filter = "true" ]; then

    # Make directory
    mkdir $work_dir/filt_VCF

    # Run filtering
    source ./vcf_filter_parallel.sh --workdir="$work_dir" \
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

if [ $impute = "true" ]; then

    # Make directory
    mkdir -p $work_dir/imp_VCF

    # Run beagle
    beagle gt=$work_dir/filt_VCF/${study_name}_production_filt.vcf.gz \
            out=$work_dir/imp_VCF/${study_name}_production_filt_imp \
            map=./SynOp_RIL906_v1.0_GBS_monotonic.map \
            nthreads=$ncores \
            window=205

    # Extract the original VCF header
    bcftools view -h $work_dir/filt_VCF/${study_name}_production_filt.vcf.gz > $work_dir/imp_VCF/original_header.vcf

    # Extract the body of the imputed VCF (no header)
    bcftools view --no-header $work_dir/imp_VCF/${study_name}_production_filt_imp.vcf.gz > $work_dir/imp_VCF/imputed_body.vcf

    # Combine the original header with the imputed body and compress
    cat $work_dir/imp_VCF/original_header.vcf $work_dir/imp_VCF/imputed_body.vcf | bgzip > $work_dir/imp_VCF/${study_name}_production_filt_imp_reheadered.vcf.gz

    # Index the final VCF using bcftools
    bcftools index --tbi $work_dir/imp_VCF/${study_name}_production_filt_imp_reheadered.vcf.gz

    # Run final filtering
    source ./vcf_filter_parallel.sh --workdir="$work_dir" \
                                    --vcfin=$work_dir/imp_VCF/${study_name}_production_filt_imp_reheadered.vcf.gz \
                                    --vcfout=$work_dir/imp_VCF/${study_name}_production_final.vcf.gz \
                                    --taxamiss=$taxamiss \
                                    --taxahet=$taxahet \
                                    --maxdep=$maxdep \
                                    --snpmiss=$snpmiss \
                                    --snphet=$snphet \
                                    --maf=$maf \
                                    --removechr=$removechr \
                                    --ncores=$ncores

fi


exit 0;
