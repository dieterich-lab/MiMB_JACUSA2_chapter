#!/bin/bash

#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=150G
#SBATCH -c 10
#SBATCH --job-name="JACUSA_test"
#SBATCH --output=Snakejob2.txt
#SBATCH --mail-user=amina.lemsara@uni-heidelberg.de
#SBATCH --partition=general

source activate snakemake
module load java
module load samtools
module load picard-tools
module load bedtools
module load R
srun snakemake --cores --unlock

srun snakemake --cores all predict_modification --rerun-incomplete
# srun snakemake --cores all visualize_pattern --rerun-incomplete
# snakemake --dag predict_modification | dot -Tsvg > dag.svg
