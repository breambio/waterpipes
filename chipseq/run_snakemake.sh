#!/bin/bash
#SBATCH -c 64
#SBATCH --mem 720GB
#SBATCH -p long,big-mem,normal,express
conda activate chipseq

snakemake --profile profile/
