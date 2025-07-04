#!/bin/bash
#SBATCH --job-name=flyphonedb2
#SBATCH --output=flyphonedb2.out
#SBATCH --error=flyphonedb2.err
#SBATCH --nodes=1
#SBATCH --time=72:00:00
#SBATCH --ntasks=24
#SBATCH --mem=120G

source /data/ebaird/miniconda3/etc/profile.d/conda.sh

conda activate R_process7

Rscript --vanilla FlyphoneDB2.r

conda deactivate