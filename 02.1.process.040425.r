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

#library(colorRamps, lib.loc = "/data/ebaird/miniconda3/envs/R_process7/lib/R/library")
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
seu <- readRDS(paste0(mainDir,"/data/seu030425.rds"))

#########################################
## 8. QC plots
#########################################
gs <- lapply(c( "nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "dblsClus", "GFP", "percent.mt", "percribo", "FBgn0260400", "FBgn0026562"), function(i) FeaturePlot(seu, features=i, order=T, ncol=3, max.cutoff="q99") & scale_colour_gradientn(colours = cols))
gs[[length(gs)+1]] <- DimPlot(seu, group.by="Phase")
gs[[length(gs)+1]] <- DimPlot(seu, group.by="sample")
gs[[length(gs)+1]] <- DimPlot(seu, group.by="SCT_snn_res.1")
gs[[length(gs)+1]] <- DimPlot(seu, group.by="scDblFinder.class")
lay <- rbind(c(1,2,3, 4),
             c(5,6, 7, 8),
             c(9, 10, 11, 11),
             c(12, 13, 14, 13))
ggexport(marrangeGrob(grobs = gs, top=NULL, layout_matrix=lay), filename = paste0(figDir, "umap.QC.pdf"), width=18, height=18)

#########################################
## 9. Remove doublets
#########################################
seu <- subset(seu, cells=colnames(seu)[seu$scDblFinder.class != "doublet"])

#########################################
## 10. Regress phase
#########################################
#install.packages("RhpcBLASctl")
library(RhpcBLASctl)
blas_set_num_threads(32)

seu <- SCTransform(seu, verbose = FALSE, vars.to.regress=c("S.Score", "G2M.Score"), return.only.var.genes = FALSE)#, conserve.memory = TRUE, vst.flavor = "v1")
#

load("/data/ebaird/scRNAseq/R_process_IRB/seu_regPhase.RData")

seu <- RunPCA(seu, verbose = FALSE)
pdf(paste0(figDir, "elbow_reg.pdf"))
ElbowPlot(seu, ndims=30)
dev.off()
#
seu <- RunUMAP(seu, dims = 1:11)
seu <- FindNeighbors(seu, dims = 1:11, verbose = FALSE)
seu <- FindClusters(seu, verbose = FALSE, resolution=1)
#
### Plot UMAPs
gs <- lapply(c( "nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "dblsClus", "GFP", "percent.mt", "percribo", "FBgn0260400", "FBgn0026562"), function(i) FeaturePlot(seu, features=i, order=T, ncol=3, max.cutoff="q99") & scale_colour_gradientn(colours = cols))
gs[[length(gs)+1]] <- DimPlot(seu, group.by="Phase")
gs[[length(gs)+1]] <- DimPlot(seu, group.by="sample")
gs[[length(gs)+1]] <- DimPlot(seu, group.by="seurat_clusters")
gs[[length(gs)+1]] <- DimPlot(seu, group.by="scDblFinder.class")
lay <- rbind(c(1,2,3, 4),
             c(5,6, 7, 8),
             c(9, 10, 11, 11),
             c(12, 13, 14, 14))
ggexport(marrangeGrob(grobs = gs, top=NULL, layout_matrix=lay), filename = paste0(figDir, "umap.QC_regPhase.pdf"), width=18, height=18)
#
pdf(paste0(figDir, "umap_sample_regPhase.pdf"), width=12, height=6)
DimPlot(seu, group.by="sample", split.by="sample")
dev.off()
#
save(seu, file=paste0(repDir, "seu_regPhase.RData"))
#
#########################################
## 10. FindMarkers
#########################################
DefaultAssay(seu)<-'SCT'
Idents(seu) <- paste0('SCT_snn_res.1')
all.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.25)
all.markers$gene <- rownames(all.markers)
all.markers$genesymbol <- mapIds(org.Dm.eg.db, keys=all.markers$gene, column="SYMBOL", keytype="FLYBASE", multiVals="first")
write.xlsx(all.markers,file=paste0(tabDir,'allMarkers.xlsx'), rowNames=TRUE)
write.csv(all.markers,file=paste0(tabDir,'allMarkers.csv'))
all.markers %>%
        group_by(cluster) %>%
        slice_max(n = 10, order_by = avg_log2FC) -> top10
write.xlsx(top10,file=paste0(tabDir,'top10Markers.xlsx'), rowNames=TRUE)
write.csv(top10,file=paste0(tabDir,'top10Markers.csv'))

#
DefaultAssay(seu)<-'SCT'
pdf(paste0(figDir,'top10markers.heatmap.pdf'),width=10,height=15)
print(DoHeatmap(seu, features = top10$gene) + NoLegend())
dev.off()
#
y<-seu
DefaultAssay(y) <- "SCT"
all.markers %>%
        group_by(cluster) %>%
        slice_max(n = 1, order_by = avg_log2FC) -> top1
pdf(paste0(figDir,'top1marker.umap.pdf'),width=20,height=15)
print(FeaturePlot(y, features = top1$gene, label=TRUE) & scale_colour_gradientn(colours =cols))
dev.off()


all.markers %>%
        group_by(cluster) %>%
        slice_max(n = 10, order_by = avg_log2FC)



geno <- c("gal", "flp", "gal", "flp")
names(geno) <- c("2196", "2197", "2198", "2199")
geno <- geno[seu$sample]
names(geno) <- colnames(seu)
seu$genotype <- geno
timepoint <- c("12d", "10d", "10d", "12d")
names(timepoint) <- c("2196", "2197", "2198", "2199")
timepoint <- timepoint[seu$sample]
names(timepoint) <- colnames(seu)
seu$timepoint <- timepoint
seu$condition <- paste0(seu$genotype, "_", seu$timepoint)

seu$condition <- paste0(geno[seu$sample], "_", timepoint[seu$sample])
seu$genotype <- geno[seu$sample]
seu$timepoint <- timepoint[seu$sample]

# cell type proportion changes

library(reshape2)
ta <- melt(t(table(seu$condition)))
colnames(ta) <- c("cluster", "condition", "ncells")
g <- ggplot(aes(condition, ncells), data=ta) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1))
pdf(paste0(figDir, "table.condition.pdf"), width=10, height=5)
print(g + theme_bw())
dev.off()
write.xlsx(ta,file=paste0(tabDir, "table.condition.xlsx"))    

ta <- melt(t(table(seu$condition,seu$SCT_snn_res.1)))
colnames(ta) <- c("condition", "cluster", "ncells")
g <- ggplot(aes(condition, ncells, fill=cluster), data=ta) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1))

g1 <- ggplot(aes(cluster, ncells, fill=cluster), data=ta) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1)) + facet_wrap(~condition, ncol=4) + scale_fill_manual(values = mycols17) + theme_bw() 
g2 <- ggplot(aes(condition, ncells, fill=condition), data=ta) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1)) + facet_wrap(~cluster, ncol=16) + scale_fill_manual(values = mycols17) + theme_bw() 


pdf(paste0(figDir, "barplot_cluster_cond.pdf"), width=10, height=5)
print(g + theme_bw())
print(g1)
print(g2)
dev.off()
write.xlsx(ta,file=paste0(tabDir, "barplot_celltype_cond.xlsx"))

#ta <- table(x$condition,x$SCT_snn_res.1)
ta <- table(seu$SCT_snn_res.1,seu$condition)
ta <- t(t(ta)/colSums(ta))
ta <- melt(t(ta))
colnames(ta) <- c("condition", "cluster", "ncells")
ta$cluster<-as.factor(ta$cluster)
g <- ggplot(aes(condition, ncells, fill=cluster), data=ta) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1) + scale_fill_manual(values = mycols17))
pdf(paste0(figDir, "barplot_cluster_cond_perc.pdf"), width=10, height=5)
print(g + theme_bw())
dev.off()
write.xlsx(ta,file= paste0(tabDir,'barplot_cluster_cond_perc.xlsx'))

# cell type proportions
ta <- table(x$SCT_snn_res.1,x$condition)
ta <- t(t(ta)/colSums(ta))
ta <- melt(t(ta))
colnames(ta) <- c("condition", "cluster", "ncells")
#ta$geno<- substr(ta$condition,6,7)
#ta$cluster<-as.factor(ta$cluster)

pdf(paste0(figDir, "cell_type_composition.pdf"), width=15, height=12)
g <- ggplot(ta,aes(geno, ncells,color=geno))+geom_point()+ geom_smooth(method='lm', fullrange=TRUE,aes(color=geno,fill=geno)) + facet_wrap(~cluster,scales = "free")+ theme_bw()
print(g)
g <- ggplot(newta,aes(geno, ncell.boxcox,color=geno))+geom_point()+ geom_smooth(method='lm', fullrange=TRUE,aes(color=geno,fill=geno)) + facet_wrap(~cluster,scales = "free")+ theme_bw() 
print(g)
dev.off()


# 10.2 cell cycle proportions

ta <- melt(t(table(seu$Phase)))
colnames(ta) <- c("cluster", "Phase", "ncells")
g <- ggplot(aes(Phase, ncells), data=ta) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1))
pdf(paste0(figDir, "table.Phase.pdf"), width=10, height=5)
print(g + theme_bw())
dev.off()
write.xlsx(ta,file=paste0(tabDir, "table.Phase.xlsx"))    

ta <- melt(t(table(seu$condition,seu$Phase)))
colnames(ta) <- c("condition", "Phase", "ncells")
g <- ggplot(aes(condition, ncells, fill=Phase), data=ta) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1))

g1 <- ggplot(aes(Phase, ncells, fill=Phase), data=ta) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1)) + facet_wrap(~condition, ncol=4) + scale_fill_manual(values = mycols17) + theme_bw() 
g2 <- ggplot(aes(condition, ncells, fill=condition), data=ta) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1)) + facet_wrap(~Phase, ncol=16) + scale_fill_manual(values = mycols17) + theme_bw() 


pdf(paste0(figDir, "barplot_Phase_cond.pdf"), width=10, height=5)
print(g + theme_bw())
print(g1)
print(g2)
dev.off()
write.xlsx(ta,file=paste0(tabDir, "barplot_Phase_cond.xlsx"))

#ta <- table(x$condition,x$SCT_snn_res.1)
ta <- table(seu$Phase,seu$condition)
ta <- t(t(ta)/colSums(ta))
ta <- melt(t(ta))
colnames(ta) <- c("condition", "Phase", "ncells")
ta$cluster<-as.factor(ta$cluster)
g <- ggplot(aes(condition, ncells, fill=cluster), data=ta) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1) + scale_fill_manual(values = mycols17))
pdf(paste0(figDir, "barplot_cluster_cond_perc.pdf"), width=10, height=5)
print(g + theme_bw())
dev.off()
write.xlsx(ta,file= paste0(tabDir,'barplot_cluster_cond_perc.xlsx'))

ta <- t(t(ta)/colSums(ta))
ta <- melt(t(ta))
colnames(ta) <- c("cluster", "Phase", "ncells")
g <- ggplot(aes(Phase, ncells), data=ta) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1))
pdf(paste0(figDir, "table.Phase.perc.pdf"), width=10, height=5)
print(g + theme_bw())
dev.off()
write.xlsx(ta,file=paste0(tabDir, "table.Phase.perc.xlsx")) 


g <- ggplot(aes(condition, ncells, fill=Phase), data=ta) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1))

g1 <- ggplot(aes(Phase, ncells, fill=Phase), data=ta) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1)) + facet_wrap(~condition, ncol=4) + scale_fill_manual(values = mycols17) + theme_bw() 
g2 <- ggplot(aes(condition, ncells, fill=condition), data=ta) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11, color=1)) + facet_wrap(~Phase, ncol=16) + scale_fill_manual(values = mycols17) + theme_bw() 


pdf(paste0(figDir, "barplot_Phase_cond_perc.pdf"), width=10, height=5)
print(g + theme_bw())
print(g1)
print(g2)
dev.off()
write.xlsx(ta,file=paste0(tabDir, "barplot_Phase_cond_perc.xlsx"))

# cell type proportions
ta <- table(x$SCT_snn_res.1,x$condition)
ta <- t(t(ta)/colSums(ta))
ta <- melt(t(ta))
colnames(ta) <- c("condition", "cluster", "ncells")
#ta$geno<- substr(ta$condition,6,7)
#ta$cluster<-as.factor(ta$cluster)

pdf(paste0(figDir, "cell_type_composition.pdf"), width=15, height=12)
g <- ggplot(ta,aes(geno, ncells,color=geno))+geom_point()+ geom_smooth(method='lm', fullrange=TRUE,aes(color=geno,fill=geno)) + facet_wrap(~cluster,scales = "free")+ theme_bw()
print(g)
g <- ggplot(newta,aes(geno, ncell.boxcox,color=geno))+geom_point()+ geom_smooth(method='lm', fullrange=TRUE,aes(color=geno,fill=geno)) + facet_wrap(~cluster,scales = "free")+ theme_bw() 
print(g)
dev.off()


#########################################
## 11. Compute MAGIC
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

load("/data/ebaird/scRNAseq/R_process_IRB/seu.RData")
geno <- c("gal", "flp", "gal", "flp")
names(geno) <- c("2196", "2197", "2198", "2199")
geno <- geno[seu$sample]
names(geno) <- colnames(seu)
seu$genotype <- geno
timepoint <- c("12d", "10d", "10d", "12d")
names(timepoint) <- c("2196", "2197", "2198", "2199")
timepoint <- timepoint[seu$sample]
names(timepoint) <- colnames(seu)
seu$timepoint <- timepoint
seu$condition <- paste0(seu$genotype, "_", seu$timepoint)

#########################################
## 12. Read signatures from Bulk RNAseq and ATACseq
#########################################
# Bulk RNAseq candidates
tmp<-data.frame(read.table("/data/ebaird/refs/o4.lrt.FLP_T0.vs.GAL_T0.csv",sep=',',header=TRUE,row.names=1))
tmp$gene<-sapply(strsplit(rownames(tmp),"_"), `[`, 1)
bulkrnaseq<-list()
tmp<-tmp[tmp$FDR<0.05,]
bulkrnaseq$FLPvsGAL_T0_UP.top100<-tmp[order(tmp$logFC),'gene'][0:100]
bulkrnaseq$FLPvsGAL_T0_DW.top100<-tmp[order(-tmp$logFC),'gene'][0:100]
save(bulkrnaseq,file=paste0(repDir,'RNAseq_candidates_20250404.RData'))

# Bulk ATACseq candidates
tmp<-read.xlsx('/data/ebaird/refs/4.Lidia.T0.GAL.vs.T0.FLP.ATAC.master.xlsx')

bulkatacseq<-list()
tmp<-tmp[tmp$FDR<0.05,]
tmp<- tmp[tmp$Nearest.PromoterID!='',]

bulkatacseq$FLPvsGAL_T0_OPEN.top100<-tmp[order(tmp$Fold, reverse=TRUE),'Nearest.PromoterID'][0:100]
bulkatacseq$FLPvsGAL_T0_CLOSE.top100<-tmp[order(tmp$Fold),'Nearest.PromoterID'][0:100]
#
save(bulkatacseq,file=paste0(repDir, 'ATAC_candidates_20250404.RData')) 

csl <- c(bulkatacseq,bulkrnaseq)

#########################################
## 13. Compute signature scores
#########################################
seu <- AddModuleScore(object = seu,features = csl, name = paste0(names(csl),"_ams"),assay='MAGIC_SCT')
colnames(seu@meta.data)[which(regexpr("_ams",colnames(seu@meta.data))>0)] <- names(csl)
seu[['MAGIC_SCT_AddModuleScore']]<-CreateAssayObject(t(seu@meta.data[,names(csl)]))

goi<-c( "FBgn0260400", "FBgn0026562")
selected<-c(names(bulkrnaseq),names(bulkatacseq))
goi_rna<-unique(unlist(bulkrnaseq))
goi_rna<- goi_rna[goi_rna %in% rownames(seu@assays$MAGIC_SCT)]
goi_atac<-unique(unlist(bulkatacseq))
goi_atac<- goi_atac[goi_atac %in% rownames(seu@assays$MAGIC_SCT)]
goi<-c(unique(goi_rna,goi_atac))

for (g in goi){
    seu[[g]]<- as.numeric(seu@assays$MAGIC_SCT[g,])
}
#
dir.create(paste0(figDir,'boxplot.MAGIC_SCT_AddModuleScore'))
dir.create(paste0(figDir,'violin.MAGIC_SCT_AddModuleScore'))

# RNAseq
tmp <- names(bulkrnaseq)
DefaultAssay(seu) <- "MAGIC_SCT_AddModuleScore"
g <- lapply(tmp, function(i) {
        cat(i, "\t")
        FeaturePlot(seu, reduction='umap',features=gsub('_','-',i), pt.size=0.5, combine=TRUE, order=TRUE,split.by='condition', label=TRUE) & scale_colour_gradientn(colours =cols)
})
g[[length(g)+1]] <- DimPlot(seu, reduction='umap', pt.size=0.5, group.by="SCT_snn_res.1",split.by='condition')
ggexport(ggarrange(plotlist = g, nrow = 10, ncol = 1), filename = paste0(figDir, "umap.bulkRNAseq.png"), width=1600, height=4000)

# ATACseq
tmp <- names(bulkatacseq)
DefaultAssay(seu) <- "MAGIC_SCT_AddModuleScore"
g <- lapply(tmp, function(i) {
        cat(i, "\t")
        FeaturePlot(seu, reduction='umap',features=gsub('_','-',i), pt.size=0.5, combine=TRUE, order=TRUE,split.by='sample') & scale_colour_gradientn(colours =cols)
})
g[[length(g)+1]] <- DimPlot(seu, reduction='umap', pt.size=0.5, group.by="SCT_snn_res.1",split.by='sample')
ggexport(ggarrange(plotlist = g, nrow = 10, ncol = 1), filename = paste0(figDir, "umap.bulkATACseq.png"), width=1600, height=4000)

# violin plots per cluster and per current annotation, split by sample
pdf(paste0(figDir, "violin.MAGIC_SCT_AddModuleScore.cell_types.pdf"), width=7.5, height=5)
for(i in names(bulkrnaseq)){
    g1 <- ggplot(seu@meta.data, aes_string(x='SCT_snn_res.1', y = i)) + geom_violin() + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    print(g1)
}
dev.off()


#group.cols<-list(
#    2198='#756bb1',
#    2197='#bcbddc',
#    2196='#de2d26',
#    2199='#fc9272')

for(i in names(bulkrnaseq)){
    df <- data.frame(expression = as.numeric(seu@assays$MAGIC_SCT_AddModuleScore@data[gsub('_','-',i),]), condition = seu$condition,cluster=seu$SCT_snn_res.1)
    png(paste0(figDir, "boxplot.MAGIC_SCT_AddModuleScore/",i,".png"), width=2000/1.75, height=3000/1.75)
    g1 <- ggplot(df, aes(x=condition, y = expression))+ geom_violin() +geom_jitter()  + theme_bw() + ggtitle(i)+facet_wrap(~cluster,ncol=4)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+ scale_fill_manual(values = group.cols)
    print(g1)
    dev.off()
}
#Also genes of interest
DefaultAssay(seu)<-'MAGIC_SCT'
for (j in names(bulkrnaseq)){
    tmp<-bulkrnaseq[[j]]
    tmp<-tmp[tmp %in% rownames(seu)]
    g <- lapply(tmp, function(i) {
        cat(i, "\t")
        FeaturePlot(seu, reduction='umap',features=gsub('_','-',i), pt.size=0.5, combine=TRUE, order=TRUE) & scale_colour_gradientn(colours =cols)
    })
    ggexport(ggarrange(plotlist = g, nrow = ceiling(length(tmp)/4), ncol = 4), filename = paste0(figDir, "umap.",j,".png"), width=1800, height=ceiling((length(tmp)+1)/4)*400)
}
#