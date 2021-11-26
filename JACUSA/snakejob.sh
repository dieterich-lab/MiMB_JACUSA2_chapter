#!/bin/bash

#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH -c 10
#SBATCH --job-name="JACUSA_HEK293T"
#SBATCH --output=JACUSA_HEK293T.txt
#SBATCH --mail-user=amina.lemsara@uni-heidelberg.de
#SBATCH --partition=general

source activate snakemake
module load java
module load samtools
module load picard-tools
module load bedtools
module load R
srun snakemake --cores --unlock

srun snakemake --cores all analysis_aggregate
