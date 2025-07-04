#!/bin/bash
#SBATCH --job-name=scVelo
#SBATCH --output=scVelo_job.out
#SBATCH --error=scVelo_job.err
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=128G

source /data/ebaird/miniconda3/etc/profile.d/conda.sh
conda activate scvelo1

# Read the output directory from the text file
MAINDIR="/data/ebaird/scRNAseq/SCENTINELsep24/"

OUTPUT_DIR="${MAINDIR}/scvelo/neuronal_subset"

# Set arguments for AnnData file and input directory
MERGED_ANNDATA_FILE="/data/ebaird/scRNAseqreports/res/Gal10d_Gal12d_Flp10d_Flp12d_070525/anndata_neuronal_subset.h5ad"
INPUT_DIR="/data/ebaird/scRNAseq/SCENTINELsep24/cellranger"
SAMPLE_NAMES=("flp_10d" "flp_12d" "gal_10d" "gal_12d")

# Construct LOOM_FILES argument dynamically
LOOM_FILES=$(for SAMPLE in "${SAMPLE_NAMES[@]}"; do echo -n "$INPUT_DIR/$SAMPLE/velocyto/$SAMPLE.loom,"; done)
LOOM_FILES=${LOOM_FILES%,} # Remove trailing comma

mkdir -p "$OUTPUT_DIR"

# Run the Python script with the arguments
python /data/ebaird/scRNAseq/SCENTINELsep24/code/scVelo.py --output_dir "$OUTPUT_DIR" --merged_anndata_file "$MERGED_ANNDATA_FILE" --loom_files "$LOOM_FILES" --sample_names "${SAMPLE_NAMES[@]}"

conda deactivate