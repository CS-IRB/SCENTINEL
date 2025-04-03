##################################################################################################
### Project Name:  ebaird_070225
### Project Type:  scell
### Species:  fly
### Description:  Data generated at IRB from brain tumors in fly
### Keywords: 
### Analyst:  ebaird
### Project Date: 070225
##################################################################################################


#########################################
## 0. Libraries and directories
#########################################

library(colorRamps, lib.loc = "/data/ebaird/miniconda3/envs/R_process7/lib/R/library")
library(Matrix)
library(parallel)
library(Seurat)
library(gridExtra)
library(ggplot2)
library(Rmagic)
library(future)
library(ggpubr)
library(viridis)
library(openxlsx)
library(viridis)
library(scran)
library(gghighlight)
library(dplyr)
library(org.Dm.eg.db)
library(ggExtra)
library(scDblFinder)
options(device=pdf)

options(future.globals.maxSize = 64912896000)

mainDir <- "/data/ebaird/scRNAseq"
dataDir <- paste0(mainDir, "data/")
repDir <- paste0(mainDir, "reports/res/02.1.process/")
figDir <- paste0(repDir, "figs/")
tabDir <- paste0(repDir, "tables/")

dir.create(repDir, recursive = TRUE, showWarnings = FALSE)
dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabDir, recursive = TRUE, showWarnings = FALSE)


mycols <- c(1, '#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506')
mycols11 <- c(1, '#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506', "purple", "violet", "gray")
mycols13 <- c(1, '#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506', "purple", "violet", "gray", "blue", "green")
mycols17 <- c(1, '#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506', "purple", "violet", "gray", "blue", "green", rainbow(4))

mycols20 <- c("yellow", '#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506', "purple", "violet", "chartreuse", "blue", "green", rainbow(4), "darkslategray3", "darksalmon", "darkorchid4", "cyan")

corner <- function(x) x[1:5,1:5]
cols <- c(colorRamps::matlab.like2(20)[1:18], "deeppink2", "deeppink3", "deeppink4")

getdensity <- function(x, y, ...) {
      dens <- MASS::kde2d(x, y, ...)
      ix <- findInterval(x, dens$x)
      iy <- findInterval(y, dens$y)
      ii <- cbind(ix, iy)
      return(dens$z[ii])
}

library(future)
plan("multicore", workers = 1)

#########################################
## 1. Import data
#########################################

samDir <- "/data/ebaird/scRNAseq/SCENTINELsep24/cellranger_res/"
sams <- dir(samDir)
samDir <- paste0(samDir, sams, "/outs/filtered_feature_bc_matrix/")
names(samDir) <- sams

#### Read matrix format

barcode.path <- paste0(samDir, "barcodes.tsv.gz")
features.path <- paste0(samDir, "features.tsv.gz")
matrix.path <- paste0(samDir, "matrix.mtx.gz")
file.exists(barcode.path)

mats <- lapply(1:length(barcode.path), function(i){
    mat <- readMM(file = matrix.path[i])
    gene.names = read.delim(features.path[i],
        header = FALSE,
        stringsAsFactors = FALSE)
    barcode.names = read.delim(barcode.path[i],
        header = FALSE,
        stringsAsFactors = FALSE)
    colnames(mat) = barcode.names$V1
    rownames(mat) = gene.names$V1
    mat
})
names(mats) <- names(samDir)

#########################################
## 2. Filtering and QC
#########################################

### Check filters

filters <- list(6000, 400, 500)
names(filters) <- c("maxGene", "minGene", "minUMI")

ncs <- lapply(mats, function(mat){
    cs <- colSums(mat)
    cg <- colSums(mat > 0)
    sapply(seq(450, to=10000, by=250), function(mu) sum(cg < filters$maxGene & cg > filters$minGene & cs > mu))
})

pdf(paste0(figDir, "filters.UMIs.pdf"))
plot(seq(450, to=10000, by=250), ncs[[1]], type="b", col=1, ylab="Number of cells passing filter", xlab="Filter on total number of UMIs per cell", ylim=c(0,12000))
for(i in 2:length(mats)) 
    lines(seq(450, to=10000, by=250), ncs[[i]], type="b", col=i, ylab="Number of cells passing filter", xlab="Filter on total number of UMIs per cell")
legend("topright", fill=1:4, legend=names(ncs))
dev.off()

ngs <- lapply(mats, function(mat){
    cs <- colSums(mat)
    cg <- colSums(mat > 0)
    sapply(seq(400, to=2000, by=50), function(mu) sum(cg < filters$maxGene & cg > mu & cs > filters$minUMI))
})

pdf(paste0(figDir, "filters.GENEs.pdf"))
plot(seq(400, to=2000, by=50), ngs[[1]], type="b", col=1, ylab="Number of cells passing filter", xlab="Filter on total number of genes per cell", ylim=c(0,12000))
for(i in 2:length(mats)) 
    lines(seq(400, to=2000, by=50), ngs[[i]], type="b", col=i)
legend("topright", fill=1:4, legend=names(ngs))
dev.off()

### Check number of genes and UMIs

dfs <- lapply(names(mats), function(i){
    mat <- mats[[i]]
    data.frame(sample=i, umis=colSums(mat), ngenes=colSums(mat > 0))
})
dfs <- do.call(rbind, dfs)


png(paste0(figDir, "umis.vs.ngenes.png"), width=1500, height=1000)
ggplot(aes(umis, ngenes), data=dfs) +  geom_bin2d(bins = 1000) + facet_wrap(.~sample)
dev.off()

pdf(paste0(figDir, "dens_umis.pdf"), width=15, height=10)
ggplot(aes(umis, fill=sample), data=dfs) + geom_density(alpha=0.5) + coord_trans(x = "log10") + geom_vline(xintercept=700)
ggplot(aes(ngenes, fill=sample), data=dfs) + geom_density(alpha=0.5) + coord_trans(x = "log10") + geom_vline(xintercept=300)
dev.off()

### Use relaxed reads for first round

sink(paste0(tabDir, "num_cells_pass.txt"))
table(dfs$umis > 700 & dfs$ngenes > 300, dfs$sample)
sink()

seus <- lapply(mats, function(mat) CreateSeuratObject(counts=mat))

seus <- lapply(seus, function(seu) subset(seu, cells=Cells(seu)[seu$nCount_RNA > 5000 & seu$nFeature_RNA > 1200]))
for(i in names(seus)) seus[[i]]$sample <- i
for(i in names(seus)) seus[[i]]$orig.ident <- NULL

### Merge into a single object

seu <- merge(seus[[1]], seus[2:4])

#########################################
## 3. Check ribosomal content
#########################################

library(rtracklayer)

gf <- readGFF("/data/ebaird/scRNAseq/R_process_IRB/genes.gtf")

syms <- gf$gene_name; names(syms) <- gf$gene_id
syms <- syms[rownames(seus[[1]])]

# library(RibosomalQC)

ribos <- names(syms)[grepl("RpL|RpS", syms)]

seur <- RibosomalQC(seu = seus, features_vec = ribos, platform = '10X')

pdf(paste0(figDir, 'histogram_RPA.pdf'))
for(i in names(seur)) print(hist(seur[[i]]))
dev.off()

pdf(paste0(figDir, 'barplots_RPA.pdf'))
for(i in names(seur)) 
print(barplot(seur[[i]], type = 3) + ggtitle(i))
dev.off()


pdf(paste0(figDir, 'quantile_plot.pdf'))
for(i in names(seur))
print(plot_quantile(seur[[i]]) + ggtitle(i))
dev.off()

for(i in names(seuf)) seuf[[i]]$percribo <- seur[[i]]$seurat$RPA

rm(seur)


#########################################
## 4. Remove RPAs and check most variable genes
#########################################
seuf <- lapply(seus, function(x) {
    x[["percribo"]] <- PercentageFeatureSet(x, features=ribos) 
x   
})

seuf <- lapply(seuf, function(x) subset(x, features=rownames(x)[!rownames(x) %in% ribos]))

seu <- merge(seuf[[1]], seuf[2:4])
seu <- JoinLayers(seu)
seu <-NormalizeData(seu)

genevar <- data.frame(modelGeneVar(seu@assays$RNA$data))
genevar$gene <- rownames(genevar)
DefaultAssay(seu) <- "RNA"

mitos <- names(syms)[grepl("^mt:", syms)]

### NO mitochondrial genes??
seu[["percent.mt"]] <- PercentageFeatureSet(seu, features=mitos)

pdf(paste0(figDir, "mean_var.pdf"), width=16)
ggplot(aes(mean,total),data=genevar)+geom_point(size=0.7)+scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')
dev.off()

#########################################
## 5. Normalize
#########################################
seu <- SCTransform(seu, verbose = FALSE, return.only.var.genes = FALSE)

seu <- RunPCA(seu, verbose = FALSE)

pdf(paste0(figDir, "elbow.pdf"))
ElbowPlot(seu, ndims=30)
dev.off()

seu <- RunUMAP(seu, dims = 1:16)
seu <- FindNeighbors(seu, dims = 1:16, verbose = FALSE)
seu <- FindClusters(seu, verbose = FALSE, resolution=1)

umapplot <- FeaturePlot(seu, features="percribo", order=T, max.cutoff="q99") & scale_colour_gradientn(colours = cols)
pdf(paste0(figDir, "umap.percribo.pdf"), width=8, height=8)
print(umapplot)
dev.off()

#########################################
## 6. Compute Cell cycle
#########################################

all_cell_cycle <- read.table("/data/ebaird/refs/Dm_cell_cycle_genes.csv", sep=",", header=T)

g2 <- all_cell_cycle[all_cell_cycle$phase=="G2/M",]$geneID
s <- all_cell_cycle[all_cell_cycle$phase=="S",]$geneID

seu <- CellCycleScoring(seu, s.features = s, g2m.features = g2, set.ident = FALSE)

#########################################
## 7. Doublets
#########################################

seuf <- SplitObject(seu, split.by="sample")

seuf <- lapply(seuf, function(x) {
    x <- SCTransform(x, verbose = FALSE, return.only.var.genes = FALSE)
    RunPCA(x, verbose = FALSE)
})

pdf(paste0(figDir, "elbows.pdf"))
for(i in names(seuf)) print(ElbowPlot(seuf[[i]], ndims=30))
dev.off()

seuf <- lapply(seuf, function(x) {
    x <- RunUMAP(x, dims = 1:16)
    x <- FindNeighbors(x, dims = 1:16, verbose = FALSE)
    FindClusters(x, verbose = FALSE, resolution=1)
})


sce <- lapply(seuf, function(x) as.SingleCellExperiment(x))

dbls <- lapply(sce, function(sc) scDblFinder(sc, clusters=NULL, BPPARAM=BiocParallel::MulticoreParam(5)))

for(i in names(dbls)) seuf[[i]]$scDblFinder.class <- dbls[[i]]$scDblFinder.class

pdf(paste0(figDir, "doublets.pdf"))
for(i in names(seuf)) print(DimPlot(seuf[[i]], group.by="scDblFinder.class"))
dev.off()

mds <- lapply(seuf, function(x) x@meta.data)
for(i in 1:4) rownames(mds[[i]]) <- paste0(rownames(mds[[i]]), "_", i)
mds <- do.call(rbind, unname(mds))

seu$scDblFinder.class <- mds[colnames(seu),]$scDblFinder.class

library(reshape2)
ta <- melt(table(seu$scDblFinder.class, seu$seurat_clusters))

ta <- data.frame(ta %>% group_by(Var2) %>% summarise(prop=value[2]/sum(value), tot=sum(value)))

seu$dblsClus <- ta$prop[match(seu$seurat_clusters, ta$Var2)]

#########################################
## 8. QC plots
#########################################


pdf(paste0(figDir, "umap.QC.pdf"), width=16, height=6)
DimPlot(seu, group.by="sample", split.by="sample") + NoLegend()
DimPlot(seu, group.by="scDblFinder.class", split.by="sample") + NoLegend()
DimPlot(seu, group.by="Phase", split.by="sample") + NoLegend()
FeaturePlot(seu, features="GFP", split.by="sample")
dev.off()

pdf(paste0(figDir, "umap.QCphasemerged.pdf"), width=8, height=8)
DimPlot(seu, group.by="Phase")
dev.off()


gs <- lapply(c( "nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "dblsClus", "GFP", "percent.mt", "FBgn0260400", "FBgn0026562"), function(i) FeaturePlot(seu, features=i, order=T, ncol=3, max.cutoff="q99") & scale_colour_gradientn(colours = cols))
gs[[length(gs)+1]] <- DimPlot(seu, group.by="Phase")
gs[[length(gs)+1]] <- DimPlot(seu, group.by="sample")
gs[[length(gs)+1]] <- DimPlot(seu, group.by="SCT_snn_res.1")
gs[[length(gs)+1]] <- DimPlot(seu, group.by="scDblFinder.class")

lay <- rbind(c(1,2,3, 4),
             c(5,6, 7, 8),
             c(9, 10, 11, 11),
             c(12, 12, 13, 13))

ggexport(marrangeGrob(grobs = gs, top=NULL, layout_matrix=lay), filename = paste0(figDir, "umap.QC.1.pdf"), width=18, height=18)

#########################################
## 9. Remove doublets
#########################################

seu <- subset(seu, cells=colnames(seu)[seu$scDblFinder.class != "doublet"])

#########################################
## 10. Regress phase
#########################################

seu <- SCTransform(seu, verbose = FALSE, vars.to.regress=c("S.Score", "G2M.Score"), return.only.var.genes = FALSE, conserve.memory = TRUE, vst.flavor = "v1")

seu <- RunPCA(seu, verbose = FALSE)
pdf(paste0(figDir, "elbow_reg.pdf"))
ElbowPlot(seu, ndims=30)
dev.off()

seu <- RunUMAP(seu, dims = 1:14)
seu <- FindNeighbors(seu, dims = 1:14, verbose = FALSE)
seu <- FindClusters(seu, verbose = FALSE, resolution=1)

### Plot UMAPs

gs <- lapply(c( "nCount_RNA", "nFeature_RNA", "percribo", "S.Score", "G2M.Score", "dblsClus", "EGFP"), function(i) FeaturePlot(seu, features=i, order=T, ncol=3, max.cutoff="q99") & scale_colour_gradientn(colours = cols))
gs[[length(gs)+1]] <- DimPlot(seu, group.by="Phase")
gs[[length(gs)+1]] <- DimPlot(seu, group.by="sample")
gs[[length(gs)+1]] <- DimPlot(seu, group.by="seurat_clusters")
gs[[length(gs)+1]] <- DimPlot(seu, group.by="scDblFinder.class")

lay <- rbind(c(1,2,3, 4),
             c(5,6, 7, 8),
             c(9, 10, 11, 11))

ggexport(marrangeGrob(grobs = gs, top=NULL, layout_matrix=lay), filename = paste0(figDir, "umap.QC_regPhase.pdf"), width=18, height=12)


pdf(paste0(figDir, "umap_sample_regPhase.pdf"), width=12, height=6)
DimPlot(seu, group.by="sample", split.by="sample")
dev.off()

#########################################
## 10. Compute MAGIC
#########################################

mag <- Rmagic::magic(seu, verbose=F)
seu@assays$MAGIC_SCT <- mag@assays$MAGIC_SCT

save(seu, file=paste0(repDir, "seu.RData"))

#########################################
## 11. Copy code to reports folder
#########################################

dir.create(paste0(repDir, "code/"))
files <- dir(); files <- files[!grepl("log", files)]
for(f in files) system(paste0("cp ", f, " ", repDir, "code/"))
