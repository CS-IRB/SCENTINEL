#!/bin/bash
#SBATCH --job-name=velocyto_run
#SBATCH --output=velocyto_run_.out
#SBATCH --error=velocyto_run.err
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --mem=240G

conda activate velocyto

# Define input directory
INPUT_DIR="/data/ebaird/scRNAseq/SCENTINELsep24/cellranger_res"
GTF_FILE="/data/ebaird/refs/genes.gtf"

velocyto run10x "$INPUT_DIR/2196" "$GTF_FILE"
velocyto run10x "$INPUT_DIR/2197" "$GTF_FILE"
velocyto run10x "$INPUT_DIR/2198" "$GTF_FILE"
velocyto run10x "$INPUT_DIR/2199" "$GTF_FILE"

conda deactivate