#!/bin/bash
#SBATCH --job-name=bam_pseudobulk
#SBATCH --output=bam_pseudobulk.out
#SBATCH --error=bam_pseudobulk.err
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=72:00:00
#SBATCH --mem=128G

OUT_DIR="/data/ebaird/scRNAseq/SCENTINELsep24/pseudobulk_output"
mkdir -p $OUT_DIR

GENOTYPES=("gal" "flp")
SAMPLES=("gal_10d" "gal_12d" "flp_10d" "flp_12d")

for geno in "${GENOTYPES[@]}"; do
  SAMPLE_FILES=$(ls /data/ebaird/scRNAseq/SCENTINELsep24/barcodes/${geno}_*.txt 2>/dev/null)
  
  if [ -z "$SAMPLE_FILES" ]; then
    echo "No samples found for $geno"
    continue
  fi

  for barcode_file in $SAMPLE_FILES; do
    samp_name=$(basename $barcode_file .txt | cut -d_ -f2-)
    sample_bam="/data/ebaird/scRNAseq/SCENTINELsep24/cellranger/${samp_name}/outs/possorted_genome_bam.bam"
    output_bam="$OUT_DIR/${samp_name}.bam"
    # Extract by CB tag using -D CB:
    samtools view -@ 4 -D CB:$barcode_file -b -o $output_bam $sample_bam
    samtools index $output_bam
  done

  ### Uncomment to Merge and index samples of same genotype or timepoint
  # merged_bam="$OUT_DIR/merged_${geno}.bam"
  # samtools merge -@ 8 $merged_bam $OUT_DIR/${geno}_*.bam
  # samtools index $merged_bam
done

### Generate BigWig files from BAM files

# for geno in "${GENOTYPES[@]}"; do
#   merged_bam="/data/ebaird/scRNAseqreports/res/Gal10d_Gal12d_Flp10d_Flp12d_070525/pseudobulk_output/merged_${geno}.bam" 
#   /data/ebaird/miniconda3/envs/deeptools/bin/bamCoverage -b $merged_bam \
#     -o $OUT_DIR/${geno}.bw \
#     --binSize 50 \
#     --normalizeUsing RPKM \
#     --numberOfProcessors 8
# done

for sample in "${SAMPLES[@]}"; do
  sample_bam="/data/ebaird/scRNAseqreports/res/Gal10d_Gal12d_Flp10d_Flp12d_070525/pseudobulk_output/*${sample}.bam" 
  /data/ebaird/miniconda3/envs/deeptools/bin/bamCoverage -b $sample_bam \
    -o $OUT_DIR/${sample}.bw \
    --binSize 50 \
    --normalizeUsing RPKM \
    --numberOfProcessors 8

done

