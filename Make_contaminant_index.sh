#!/usr/bin/env bash
#SBATCH --job-name=MakeContamIndex                  # Job name
#SBATCH --output=MakeContamIndex.%A.out             # File to which stdout will be written
#SBATCH --error=MakeContamIndex.%A.err              # File to which stderr will be written
#SBATCH --partition=2tcourt                         # Partition
#SBATCH --cpus-per-task=24                          # Number of cpus ask per task
#SBATCH --time=00-01:00                             # Runtime in DD-HH:MM
#SBATCH --mem=10G                                   # Memory for all cores in Gbytes
#SBATCH --mail-type=ALL                             # BEGIN,END,FAIL,ALL
#SBATCH --mail-user=nicolas.nesi@unicaen.fr         # Email address

# ---------------------------------
# Version v1.0
# ---------------------------------
# author: Nicolas Nesi
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date: 26/08/2026
# ---------------------------------

# ---------------------------------
# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used:$THREADS"
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
module load py_env/miniconda/25.7.0
eval "$(conda shell.bash hook)"
conda activate ptp_env
# include:
# bowtie2=2.5.5
# datasets=18.36.0
# seqkit=2.13.0
# ---------------------------------

# ---------------------------------
tail -n +2 decontamination_database.csv | tr -d '\r' > contam_list_clean.csv

grep 'NC_012920.1' contam_list_clean.csv > seq_lines.csv
grep -v 'UniVec\|NC_012920.1' contam_list_clean.csv > assembly_lines.csv

while IFS=, read -r group genome accession; do
    n=0
    until [ $n -ge 3 ]; do
        datasets download genome accession "$accession" --include genome --filename "${accession}.zip" && break
        n=$((n+1))
        echo "Retry $n for $accession"
        sleep 5
    done
    unzip -o "${accession}.zip" -d "${accession}_dir" \
        && fna_file=$(find "${accession}_dir" -name "*.fna") \
        && seqkit replace -p '^(\S+)' -r "${genome}_\$1" "$fna_file" >> genomes_renamed.fa
done < assembly_lines.csv

# accessions de séquence individuelle (efetch, ex: NC_012920.1)
while IFS=, read -r group genome accession; do
    wget -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${accession}&rettype=fasta&retmode=text" \
        | seqkit replace -p '^(\S+)' -r "${genome}_\$1" - >> genomes_renamed.fa
done < seq_lines.csv

wget -O univec.fa https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec
seqkit replace -p 'gnl\|uv\|(\S+)' -r 'UniVec_$1' univec.fa >> genomes_renamed.fa

echo "start indexing"

bowtie2-build --threads "${THREADS}" genomes_renamed.fa contam_bt2_index
# ---------------------------------