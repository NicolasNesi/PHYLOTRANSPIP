#!/usr/bin/env bash
#SBATCH --job-name=GetSraList                       # Job name
#SBATCH --output=GetSraList.%A.out                  # File to which stdout will be written
#SBATCH --error=GetSraList.%A.err                   # File to which stderr will be written
#SBATCH --partition=2tcourt                         # Partition
#SBATCH --cpus-per-task=1                           # Number of cpus ask per task
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
# Submission:
# sbatch GetSraList.sh -o $ORGANISM -s $STRATEGY -db $DATABASE
# ---------------------------------

# ---------------------------------
# Variables and Paths
# ---------------------------------
# safety for failing steps
set -euo pipefail

# Positional parameters
GROUP="ViroCaen"
PROJECT="PHYLOTRANSPIP"
ORGANISM=""
STRATEGY=""
DATABASE=""

usage() {
    echo "Usage: $0 -o ORGANISM -s STRATEGY -db DATABASE" >&2
	echo "   -o, --organism      Mandatory. e.g. Phyllostomidae" >&2
	echo "   -s, --strategy      Mandatory. e.g. RNA-Seq" >&2
	echo "   -db, --database     Mandatory. Databse selected, either sra, ena or geo" >&2
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
    -o|--organism) ORGANISM="$2"; shift 2 ;;
    -s|--strategy) STRATEGY="$2"; shift 2 ;;
    -db|--database) DATABASE="$2"; shift 2 ;;
    *) echo "Invalid option: $1" >&2; usage;;
    esac
done

if [[ -z "$ORGANISM" || -z "$STRATEGY" || -z "$DATABASE" ]]; then
    echo "Error: -o, -s, and -db are all mandatory" >&2
    usage
    exit 1
fi

# Paths
FOLDER="/dlocal/home/2019013/Data/$GROUP/$PROJECT"

RUNDATE="$(date '+%Y-%m-%d')"

# Echo variables
echo "Scientific name of the sample organism: $ORGANISM"
echo "Library preperation strategy: $STRATEGY"
echo "Database selected: $DATABASE"
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
module load py_env/miniconda/25.7.0
eval "$(conda shell.bash hook)"
conda activate ptp_env
# Include:
# pysradb=2.5.1
# ---------------------------------

# ---------------------------------
# Get list of SRA of your species of interest
# ---------------------------------
# https://www.ncbi.nlm.nih.gov/sra
# Sequence Read Archive (SRA) (NCBI) / European Nucleotide Archive (ENA) / Gene Expresion Omnibus (GEO) (NCBI)

# study accession (SRP, ERP, DRP)
# experiment accession (SRX, ERX, DRX)
# sample accession (SRS, ERS, DRS)
# run accession (SRR, ERR, DRR)

pysradb search \
--db "${DATABASE}" \
--organism "${ORGANISM}" \
--strategy "${STRATEGY}" \
--layout PAIRED \
--source TRANSCRIPTOMIC \
--selection cDNA \
--max 20000 \
--stats \
--saveto "${FOLDER}"/Scripts/"${ORGANISM}_${STRATEGY}_${DATABASE}_${RUNDATE}.tsv"

echo "SRA list saved to ${FOLDER}/Scripts/${ORGANISM}_${STRATEGY}_${DATABASE}_${RUNDATE}.tsv"

echo "Select your data of interest and create a list_sra.csv input file with the format: Genus_species_locality_code,SRRXXXXX"
# ---------------------------------