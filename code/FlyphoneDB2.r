library(FlyPhone)

RunFlyPhone(
    knowledgebase_version = "Version 2 High/Moderate",
    seuratObject = "/data/ebaird/scRNAseqreports/res/Gal10d_Gal12d_Flp10d_Flp12d_070525/seurat_with_regulons.rds",
    base_output_dir = "/data/ebaird/scRNAseq/FlyPhone_output/Gal10d_Gal12d_Flp10d_Flp12d_db2",
    control_name = "gal",
    mutant_name = "flp"
)