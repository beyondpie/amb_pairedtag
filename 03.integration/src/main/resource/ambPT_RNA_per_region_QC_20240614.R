## Generation of individual brain region seurat objects for Paired-Tag RNA modality, and QC check
## Zhaoning Wang ver.20240614


library(trqwe)
#library(BPCells)
library(EnsDb.Mmusculus.v79)
library(tibble)
library(Seurat)
library(Signac)
library(ggplot2)
library(ggforce)
library(ggrepel)
library(ggridges)
library(ggseqlogo) 
library(cowplot)
library(dplyr)
library(RColorBrewer)
library(S4Vectors)
library(sctransform)
library(pheatmap)
library(harmony)


setwd("~/Desktop/mouse_brain_paired_tag/2024_new_objects/")


brain.exp1.rna <- mcreadRDS(file = "./brain.exp1.rna.rds")
brain.exp2.rna <- mcreadRDS(file = "./brain.exp2.rna.rds")
brain.exp3.rna <- mcreadRDS(file = "./brain.exp3.rna.rds")
brain.exp4.rna <- mcreadRDS(file = "./brain.exp4.rna.rds")
brain.exp5.rna <- mcreadRDS(file = "./brain.exp5.rna.rds")
brain.exp6.rna <- mcreadRDS(file = "./brain.exp6.rna.rds")
brain.exp7.rna <- mcreadRDS(file = "./brain.exp7.rna.rds")
brain.exp8.rna <- mcreadRDS(file = "./brain.exp8.rna.rds")
brain.exp9.rna <- mcreadRDS(file = "./brain.exp9.rna.rds")
brain.exp10.rna <- mcreadRDS(file = "./brain.exp10.rna.rds")
brain.exp11.rna <- mcreadRDS(file = "./brain.exp11.rna.rds")
brain.exp12.rna <- mcreadRDS(file = "./brain.exp12.rna.rds")
brain.exp13.rna <- mcreadRDS(file = "./brain.exp13.rna.rds")
brain.exp14.rna <- mcreadRDS(file = "./brain.exp14.rna.rds")



brain <- merge(brain.exp1.rna, y=c(brain.exp2.rna, 
                                   brain.exp3.rna, 
                                   brain.exp4.rna, 
                                   brain.exp5.rna, 
                                   brain.exp6.rna, 
                                   brain.exp7.rna, 
                                   brain.exp8.rna, 
                                   brain.exp9.rna, 
                                   brain.exp10.rna, 
                                   brain.exp11.rna, 
                                   brain.exp12.rna, 
                                   brain.exp13.rna, 
                                   brain.exp14.rna), 
               project = "brain")


table(brain@meta.data$brainregion, brain@meta.data$modality)

brain

#brain <- JoinLayers(brain)  # Error: vector

mcsaveRDS(brain, file = "./brain.all.rna.rds")

rm(brain.exp1.rna,brain.exp2.rna,brain.exp3.rna,brain.exp4.rna,brain.exp5.rna,brain.exp6.rna,brain.exp7.rna,brain.exp8.rna,brain.exp9.rna,brain.exp10.rna,brain.exp11.rna,brain.exp12.rna,brain.exp13.rna,brain.exp14.rna)

gc()


#brain <- readRDS(file = "./brain.all.rna.rds")

brain <- mcreadRDS(file = "./brain.all.rna.rds")



table(brain@meta.data$brainregion, brain@meta.data$modality)

Idents(brain)<-"brainregion"

AMY <- subset(brain, idents = c("AMY"))
CPU <- subset(brain, idents = c("CPU"))
ERC <- subset(brain, idents = c("ERC"))
HCa <- subset(brain, idents = c("HCa"))
HCp <- subset(brain, idents = c("HCp"))
HYP <- subset(brain, idents = c("HYP"))
NAC <- subset(brain, idents = c("NAC"))
PFC <- subset(brain, idents = c("PFC"))
VTA_SnR <- subset(brain, idents = c("VTA_SnR"))


rm(brain)

gc()

####################
AMY <- JoinLayers(AMY)
CPU <- JoinLayers(CPU)
ERC <- JoinLayers(ERC)
HCa <- JoinLayers(HCa)
HCp <- JoinLayers(HCp)
HYP <- JoinLayers(HYP)
NAC <- JoinLayers(NAC)
PFC <- JoinLayers(PFC)
VTA_SnR <- JoinLayers(VTA_SnR)


####################

AMY <- NormalizeData(AMY, normalization.method = "LogNormalize", scale.factor = 10000)


AMY <- FindVariableFeatures(AMY, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(AMY), 10)
VariableFeaturePlot(AMY)

AMY <- ScaleData(AMY, vars.to.regress = "percent.mt")
AMY <- RunPCA(AMY, features = VariableFeatures(object = AMY)) 
ElbowPlot(AMY, ndims = 50)

DimPlot(AMY, reduction = "pca")

DimHeatmap(AMY, dims = 1:15, cells = 500, balanced = TRUE)

# AMY <- JackStraw(AMY, num.replicate = 100)
# AMY<- ScoreJackStraw(AMY, dims = 1:20)
# JackStrawPlot(AMY, dims = 1:20)
# 
AMY <- FindNeighbors(AMY,reduction="pca", dims= 1:25, k.param=25)
AMY <- FindClusters(AMY, resolution = 0.5, algorithm = 2)

AMY <- RunUMAP(AMY, dims = 1:15, metric = "euclidean", min.dist = 0.05, local.connectivity=10, spread=1, seed.use=131, umap.method = "uwot", n.neighbors = 25, uwot.sgd=TRUE, verbose=TRUE)
p <- DimPlot(AMY, pt.size=0.10)

pdf("./AMY.UMAP.pdf")
print(p)
dev.off()

SaveSeuratRds(object = AMY,file = "./AMY.rna.rds")

rm(AMY)

####################

CPU <- NormalizeData(CPU, normalization.method = "LogNormalize", scale.factor = 10000)


CPU <- FindVariableFeatures(CPU, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(CPU), 10)
VariableFeaturePlot(CPU)

CPU <- ScaleData(CPU, vars.to.regress = "percent.mt")
CPU <- RunPCA(CPU, features = VariableFeatures(object = CPU))  

ElbowPlot(CPU, ndims = 50)

DimPlot(CPU, reduction = "pca")

DimHeatmap(CPU, dims = 1:15, cells = 500, balanced = TRUE)

# CPU <- JackStraw(CPU, num.replicate = 100)
# CPU<- ScoreJackStraw(CPU, dims = 1:20)
# JackStrawPlot(CPU, dims = 1:20)
# 
CPU <- FindNeighbors(CPU,reduction="pca", dims= 1:25, k.param=25)
CPU <- FindClusters(CPU, resolution = 0.5, algorithm = 2)

CPU <- RunUMAP(CPU, dims = 1:15, metric = "euclidean", min.dist = 0.05, local.connectivity=10, spread=1, seed.use=131, umap.method = "uwot", n.neighbors = 25, uwot.sgd=TRUE, verbose=TRUE)
p<- DimPlot(CPU, pt.size=0.10)

pdf("./CPU.UMAP.pdf")
print(p)
dev.off()

SaveSeuratRds(object = CPU,file = "./CPU.rna.rds")

rm(CPU)


####################

ERC <- NormalizeData(ERC, normalization.method = "LogNormalize", scale.factor = 10000)


ERC <- FindVariableFeatures(ERC, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(ERC), 10)
VariableFeaturePlot(ERC)

ERC <- ScaleData(ERC, vars.to.regress = "percent.mt")
ERC <- RunPCA(ERC, features = VariableFeatures(object = ERC))  

ElbowPlot(ERC, ndims = 50)

DimPlot(ERC, reduction = "pca")

DimHeatmap(ERC, dims = 1:15, cells = 500, balanced = TRUE)

# ERC <- JackStraw(ERC, num.replicate = 100)
# ERC<- ScoreJackStraw(ERC, dims = 1:20)
# JackStrawPlot(ERC, dims = 1:20)
# 
ERC <- FindNeighbors(ERC,reduction="pca", dims= 1:25, k.param=25)
ERC <- FindClusters(ERC, resolution = 0.5, algorithm = 2)

ERC <- RunUMAP(ERC, dims = 1:15, metric = "euclidean", min.dist = 0.05, local.connectivity=10, spread=1, seed.use=131, umap.method = "uwot", n.neighbors = 25, uwot.sgd=TRUE, verbose=TRUE)
p <- DimPlot(ERC, pt.size=0.10)

pdf("./ERC.UMAP.pdf")
print(p)
dev.off()

SaveSeuratRds(object = ERC,file = "./ERC.rna.rds")

rm(ERC)

####################

HCa <- NormalizeData(HCa, normalization.method = "LogNormalize", scale.factor = 10000)


HCa <- FindVariableFeatures(HCa, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(HCa), 10)
VariableFeaturePlot(HCa)

HCa <- ScaleData(HCa, vars.to.regress = "percent.mt")
HCa <- RunPCA(HCa, features = VariableFeatures(object = HCa)) 

ElbowPlot(HCa, ndims = 50)

DimPlot(HCa, reduction = "pca")

DimHeatmap(HCa, dims = 1:15, cells = 500, balanced = TRUE)

# HCa <- JackStraw(HCa, num.replicate = 100)
# HCa<- ScoreJackStraw(HCa, dims = 1:20)
# JackStrawPlot(HCa, dims = 1:20)
# 
HCa <- FindNeighbors(HCa,reduction="pca", dims= 1:25, k.param=25)
HCa <- FindClusters(HCa, resolution = 0.5, algorithm = 2)

HCa <- RunUMAP(HCa, dims = 1:15, metric = "euclidean", min.dist = 0.05, local.connectivity=10, spread=1, seed.use=131, umap.method = "uwot", n.neighbors = 25, uwot.sgd=TRUE, verbose=TRUE)
p <- DimPlot(HCa, pt.size=0.10)

pdf("./HCa.UMAP.pdf")
print(p)
dev.off()

SaveSeuratRds(object = HCa,file = "./HCa.rna.rds")

rm(HCa)

####################

HCp <- NormalizeData(HCp, normalization.method = "LogNormalize", scale.factor = 10000)


HCp <- FindVariableFeatures(HCp, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(HCp), 10)
VariableFeaturePlot(HCp)

HCp <- ScaleData(HCp, vars.to.regress = "percent.mt")
HCp <- RunPCA(HCp, features = VariableFeatures(object = HCp))

ElbowPlot(HCp, ndims = 50)

DimPlot(HCp, reduction = "pca")

DimHeatmap(HCp, dims = 1:15, cells = 500, balanced = TRUE)

# HCp <- JackStraw(HCp, num.replicate = 100)
# HCp<- ScoreJackStraw(HCp, dims = 1:20)
# JackStrawPlot(HCp, dims = 1:20)
# 
HCp <- FindNeighbors(HCp,reduction="pca", dims= 1:25, k.param=25)
HCp <- FindClusters(HCp, resolution = 0.5, algorithm = 2)

HCp <- RunUMAP(HCp, dims = 1:15, metric = "euclidean", min.dist = 0.05, local.connectivity=10, spread=1, seed.use=131, umap.method = "uwot", n.neighbors = 25, uwot.sgd=TRUE, verbose=TRUE)
p <- DimPlot(HCp, pt.size=0.10)

pdf("./HCp.UMAP.pdf")
print(p)
dev.off()

SaveSeuratRds(object = HCp,file = "./HCp.rna.rds")

rm(HCp)

####################

HYP <- NormalizeData(HYP, normalization.method = "LogNormalize", scale.factor = 10000)


HYP <- FindVariableFeatures(HYP, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(HYP), 10)
VariableFeaturePlot(HYP)

HYP <- ScaleData(HYP, vars.to.regress = "percent.mt")
HYP <- RunPCA(HYP, features = VariableFeatures(object = HYP))  

ElbowPlot(HYP, ndims = 50)

DimPlot(HYP, reduction = "pca")

DimHeatmap(HYP, dims = 1:15, cells = 500, balanced = TRUE)

# HYP <- JackStraw(HYP, num.replicate = 100)
# HYP<- ScoreJackStraw(HYP, dims = 1:20)
# JackStrawPlot(HYP, dims = 1:20)
# 
HYP <- FindNeighbors(HYP,reduction="pca", dims= 1:25, k.param=25)
HYP <- FindClusters(HYP, resolution = 0.5, algorithm = 2)

HYP <- RunUMAP(HYP, dims = 1:15, metric = "euclidean", min.dist = 0.05, local.connectivity=10, spread=1, seed.use=131, umap.method = "uwot", n.neighbors = 25, uwot.sgd=TRUE, verbose=TRUE)
p <- DimPlot(HYP, pt.size=0.10)

pdf("./HYP.UMAP.pdf")
print(p)
dev.off()

SaveSeuratRds(object = HYP,file = "./HYP.rna.rds")

rm(HYP)

####################

NAC <- NormalizeData(NAC, normalization.method = "LogNormalize", scale.factor = 10000)


NAC <- FindVariableFeatures(NAC, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(NAC), 10)
VariableFeaturePlot(NAC)

NAC <- ScaleData(NAC, vars.to.regress = "percent.mt")
NAC <- RunPCA(NAC, features = VariableFeatures(object = NAC)) 
ElbowPlot(NAC, ndims = 50)

DimPlot(NAC, reduction = "pca")

DimHeatmap(NAC, dims = 1:15, cells = 500, balanced = TRUE)

# NAC <- JackStraw(NAC, num.replicate = 100)
# NAC<- ScoreJackStraw(NAC, dims = 1:20)
# JackStrawPlot(NAC, dims = 1:20)
# 
NAC <- FindNeighbors(NAC,reduction="pca", dims= 1:25, k.param=25)
NAC <- FindClusters(NAC, resolution = 0.5, algorithm = 2)

NAC <- RunUMAP(NAC, dims = 1:15, metric = "euclidean", min.dist = 0.05, local.connectivity=10, spread=1, seed.use=131, umap.method = "uwot", n.neighbors = 25, uwot.sgd=TRUE, verbose=TRUE)
p <- DimPlot(NAC, pt.size=0.10)


pdf("./NAC.UMAP.pdf")
print(p)
dev.off()

SaveSeuratRds(object = NAC,file = "./NAC.rna.rds")

rm(NAC)

####################

PFC <- NormalizeData(PFC, normalization.method = "LogNormalize", scale.factor = 10000)


PFC <- FindVariableFeatures(PFC, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(PFC), 10)
VariableFeaturePlot(PFC)

PFC <- ScaleData(PFC, vars.to.regress = "percent.mt")
PFC <- RunPCA(PFC, features = VariableFeatures(object = PFC)) 
ElbowPlot(PFC, ndims = 50)

DimPlot(PFC, reduction = "pca")

DimHeatmap(PFC, dims = 1:15, cells = 500, balanced = TRUE)

# PFC <- JackStraw(PFC, num.replicate = 100)
# PFC<- ScoreJackStraw(PFC, dims = 1:20)
# JackStrawPlot(PFC, dims = 1:20)
# 
PFC <- FindNeighbors(PFC,reduction="pca", dims= 1:25, k.param=25)
PFC <- FindClusters(PFC, resolution = 0.5, algorithm = 2)

PFC <- RunUMAP(PFC, dims = 1:15, metric = "euclidean", min.dist = 0.05, local.connectivity=10, spread=1, seed.use=131, umap.method = "uwot", n.neighbors = 25, uwot.sgd=TRUE, verbose=TRUE)
p <- DimPlot(PFC, pt.size=0.10)

pdf("./PFC.UMAP.pdf")
print(p)
dev.off()

SaveSeuratRds(object = PFC,file = "./PFC.rna.rds")

rm(PFC)

####################

VTA_SnR <- NormalizeData(VTA_SnR, normalization.method = "LogNormalize", scale.factor = 10000)


VTA_SnR <- FindVariableFeatures(VTA_SnR, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(VTA_SnR), 10)
VariableFeaturePlot(VTA_SnR)

VTA_SnR <- ScaleData(VTA_SnR, vars.to.regress = "percent.mt")
VTA_SnR <- RunPCA(VTA_SnR, features = VariableFeatures(object = VTA_SnR)) 

ElbowPlot(VTA_SnR, ndims = 50)

DimPlot(VTA_SnR, reduction = "pca")

DimHeatmap(VTA_SnR, dims = 1:15, cells = 500, balanced = TRUE)

# VTA_SnR <- JackStraw(VTA_SnR, num.replicate = 100)
# VTA_SnR<- ScoreJackStraw(VTA_SnR, dims = 1:20)
# JackStrawPlot(VTA_SnR, dims = 1:20)
# 
VTA_SnR <- FindNeighbors(VTA_SnR,reduction="pca", dims= 1:25, k.param=25)
VTA_SnR <- FindClusters(VTA_SnR, resolution = 0.5, algorithm = 2)

VTA_SnR <- RunUMAP(VTA_SnR, dims = 1:15, metric = "euclidean", min.dist = 0.05, local.connectivity=10, spread=1, seed.use=131, umap.method = "uwot", n.neighbors = 25, uwot.sgd=TRUE, verbose=TRUE)
p <- DimPlot(VTA_SnR, pt.size=0.10)

pdf("./VTA_SnR.UMAP.pdf")
print(p)
dev.off()

SaveSeuratRds(object = VTA_SnR,file = "./VTA_SnR.rna.rds")

rm(VTA_SnR)
gc()

###############################################################
setwd("~/Desktop/20240608_check_anno/")


#all.metadata <- read.csv("./20240603_pairedtag.cell.meta.all.with.neutfbyregion.csv")
all.newmetadata <- read.csv("./20240614_pairedtag.cell.meta.all.with.tfv3.csv")
# 2,722,889 total barcodes
all.metadata <- all.newmetadata
rm(all.newmetadata)

AMY <- mcreadRDS(file = "./AMY.rna.rds")

#AMY@meta.data$orig.ident
#AMY@meta.data$nCount_RNA
#AMY@meta.data$nFeature_RNA
#AMY@meta.data$percent.mt
#AMY@meta.data$brainregion
#AMY@meta.data$modality
#AMY@meta.data$oriBarcode
#AMY@meta.data$sublib
#AMY@meta.data$sex
#AMY@meta.data$rep
#AMY@meta.data$exp
## AMY@meta.data$RNA_snn_res.0.5
## AMY@meta.data$seurat_clusters


rownames(all.metadata) <- all.metadata$barcode
all.metadata["barcode"] <- NULL
all.metadata["orig.ident"] <- NULL
all.metadata["nCount_RNA"] <- NULL
all.metadata["nFeature_RNA"] <- NULL
all.metadata["percent.mt"] <- NULL
all.metadata["brainregion"] <- NULL
all.metadata["modality"] <- NULL
all.metadata["oriBarcode"] <- NULL
all.metadata["sublib"] <- NULL
all.metadata["sex"] <- NULL
all.metadata["rep"] <- NULL
all.metadata["exp"] <- NULL


all.metadata <- all.metadata[all.metadata$L5r != "LQ_tf_by_region", ]
# after this filtering, there are 2,704,133 barcodes


tmp.meta <- all.metadata[intersect(rownames(all.metadata), rownames(AMY@meta.data)),c(1:30)]

AMY
# 54287 features across 266446 samples within 1 assay 

AMY.clean <- subset(AMY, cells = rownames(tmp.meta))

AMY.clean
# 54287 features across 258803 samples within 1 assay 

AMY.clean <- AddMetaData(AMY.clean, metadata = tmp.meta)

DimPlot(AMY.clean, pt.size = 0.1, group.by = "class_id_label", raster = FALSE, label = TRUE) + NoLegend()

table(AMY.clean@meta.data$class_id_label)

DimPlot(AMY.clean, pt.size = 0.1, group.by = "subclass_id_label", raster = FALSE, label = TRUE) + NoLegend()
FeaturePlot(AMY.clean, features = c("Snap25"),cols = c("lightgrey","red2"))

FeaturePlot(AMY.clean, features = c("Slc17a7"),cols = c("lightgrey","red2"))
FeaturePlot(AMY.clean, features = c("Slc17a8"),cols = c("lightgrey","red2"))
FeaturePlot(AMY.clean, features = c("Gad1"),cols = c("lightgrey","red2"))
FeaturePlot(AMY.clean, features = c("Gad2"),cols = c("lightgrey","red2"))
FeaturePlot(AMY.clean, features = c("Atp1a2"),cols = c("lightgrey","red2"))
FeaturePlot(AMY.clean, features = c("Aqp4"),cols = c("lightgrey","red2"))
FeaturePlot(AMY.clean, features = c("Ctss"),cols = c("lightgrey","red2"))

mcsaveRDS(AMY.clean, file = "./AMY.clean.rna.rds")

#AMY.clean <- mcreadRDS(file = "./AMY.clean.rna.rds")

rm(AMY.clean, AMY, tmp.meta)

############################################################

CPU <- mcreadRDS(file = "./CPU.rna.rds")

tmp.meta <- all.metadata[intersect(rownames(all.metadata), rownames(CPU@meta.data)),c(1:30)]

CPU
# 54567 features across 360889 samples within 1 assay

CPU.clean <- subset(CPU, cells = rownames(tmp.meta))

CPU.clean
# 54567 features across 349466 samples within 1 assay 

CPU.clean <- AddMetaData(CPU.clean, metadata = tmp.meta)

DimPlot(CPU.clean, pt.size = 0.1, group.by = "class_id_label", raster = FALSE, label = TRUE) + NoLegend()

table(CPU.clean@meta.data$class_id_label)
table(CPU.clean@meta.data$subclass_id_label)

DimPlot(CPU.clean, pt.size = 0.1, group.by = "subclass_id_label", raster = FALSE, label = TRUE) + NoLegend()

FeaturePlot(CPU.clean, features = c("Snap25"),cols = c("lightgrey","red2"))

FeaturePlot(CPU.clean, features = c("Flt1"),cols = c("lightgrey","red2"))
FeaturePlot(CPU.clean, features = c("Pecam1"),cols = c("lightgrey","red2"))

FeaturePlot(CPU.clean, features = c("Slc17a7"),cols = c("lightgrey","red2"))
FeaturePlot(CPU.clean, features = c("Slc17a8"),cols = c("lightgrey","red2"))
FeaturePlot(CPU.clean, features = c("Gad1"),cols = c("lightgrey","red2"))
FeaturePlot(CPU.clean, features = c("Gad2"),cols = c("lightgrey","red2"))
FeaturePlot(CPU.clean, features = c("Atp1a2"),cols = c("lightgrey","red2"))
FeaturePlot(CPU.clean, features = c("Aqp4"),cols = c("lightgrey","red2"))
FeaturePlot(CPU.clean, features = c("Ctss"),cols = c("lightgrey","red2"))

mcsaveRDS(CPU.clean, file = "./CPU.clean.rna.rds")

#CPU.clean <- mcreadRDS(file = "./CPU.clean.rna.rds")

rm(CPU.clean, CPU, tmp.meta)


############################################################

ERC <- mcreadRDS(file = "./ERC.rna.rds")

tmp.meta <- all.metadata[intersect(rownames(all.metadata), rownames(ERC@meta.data)),c(1:30)]

ERC
# 54287 features across 339,507 samples within 1 assay

ERC.clean <- subset(ERC, cells = rownames(tmp.meta))

ERC.clean
# 54287 features across 329,349 samples within 1 assay

ERC.clean <- AddMetaData(ERC.clean, metadata = tmp.meta)

DimPlot(ERC.clean, pt.size = 0.1, group.by = "class_id_label", raster = FALSE, label = TRUE) + NoLegend()

table(ERC.clean@meta.data$class_id_label)
table(ERC.clean@meta.data$subclass_id_label)

DimPlot(ERC.clean, pt.size = 0.1, group.by = "subclass_id_label", raster = FALSE, label = TRUE) + NoLegend()
DimPlot(ERC.clean, pt.size = 0.1, group.by = "supertype_id_label", raster = FALSE, label = TRUE) + NoLegend()


FeaturePlot(ERC.clean, features = c("Snap25"),cols = c("lightgrey","red2"))

FeaturePlot(ERC.clean, features = c("Flt1"),cols = c("lightgrey","red2"))
FeaturePlot(ERC.clean, features = c("Pecam1"),cols = c("lightgrey","red2"))

FeaturePlot(ERC.clean, features = c("Slc17a7"),cols = c("lightgrey","red2"))
FeaturePlot(ERC.clean, features = c("Slc17a8"),cols = c("lightgrey","red2"))
FeaturePlot(ERC.clean, features = c("Gad1"),cols = c("lightgrey","red2"))
FeaturePlot(ERC.clean, features = c("Gad2"),cols = c("lightgrey","red2"))
FeaturePlot(ERC.clean, features = c("Atp1a2"),cols = c("lightgrey","red2"))
FeaturePlot(ERC.clean, features = c("Aqp4"),cols = c("lightgrey","red2"))
FeaturePlot(ERC.clean, features = c("Ctss"),cols = c("lightgrey","red2"))

mcsaveRDS(ERC.clean, file = "./ERC.clean.rna.rds")

#ERC.clean <- mcreadRDS(file = "./ERC.clean.rna.rds")

rm(ERC.clean, ERC, tmp.meta)

############################################################

HCa <- mcreadRDS(file = "./HCa.rna.rds")

tmp.meta <- all.metadata[intersect(rownames(all.metadata), rownames(HCa@meta.data)),c(1:30)]

HCa
# 54598 features across 431463 samples within 1 assay 

HCa.clean <- subset(HCa, cells = rownames(tmp.meta))

HCa.clean
# 54598 features across 411887 samples within 1 assay 

HCa.clean <- AddMetaData(HCa.clean, metadata = tmp.meta)

DimPlot(HCa.clean, pt.size = 0.1, group.by = "class_id_label", raster = FALSE, label = TRUE) + NoLegend()

table(HCa.clean@meta.data$class_id_label)
table(HCa.clean@meta.data$subclass_id_label)

DimPlot(HCa.clean, pt.size = 0.1, group.by = "subclass_id_label", raster = FALSE, label = TRUE) + NoLegend()
DimPlot(HCa.clean, pt.size = 0.1, group.by = "supertype_id_label", raster = FALSE, label = TRUE) + NoLegend()


FeaturePlot(HCa.clean, features = c("Snap25"),cols = c("lightgrey","red2"))

FeaturePlot(HCa.clean, features = c("Flt1"),cols = c("lightgrey","red2"))
FeaturePlot(HCa.clean, features = c("Pecam1"),cols = c("lightgrey","red2"))

FeaturePlot(HCa.clean, features = c("Slc17a7"),cols = c("lightgrey","red2"))
FeaturePlot(HCa.clean, features = c("Slc17a8"),cols = c("lightgrey","red2"))
FeaturePlot(HCa.clean, features = c("Gad1"),cols = c("lightgrey","red2"))
FeaturePlot(HCa.clean, features = c("Gad2"),cols = c("lightgrey","red2"))
FeaturePlot(HCa.clean, features = c("Atp1a2"),cols = c("lightgrey","red2"))
FeaturePlot(HCa.clean, features = c("Aqp4"),cols = c("lightgrey","red2"))
FeaturePlot(HCa.clean, features = c("Ctss"),cols = c("lightgrey","red2"))
FeaturePlot(HCa.clean, features = c("Pdgfra"),cols = c("lightgrey","red2"))
FeaturePlot(HCa.clean, features = c("Mbp"),cols = c("lightgrey","red2"))

mcsaveRDS(HCa.clean, file = "./HCa.clean.rna.rds")

#HCa.clean <- mcreadRDS(file = "./HCa.clean.rna.rds")

rm(HCa.clean, HCa, tmp.meta)


############################################################

HCp <- mcreadRDS(file = "./HCp.rna.rds")

tmp.meta <- all.metadata[intersect(rownames(all.metadata), rownames(HCp@meta.data)),c(1:30)]

HCp
# 54581 features across 319140 samples within 1 assay 

HCp.clean <- subset(HCp, cells = rownames(tmp.meta))

HCp.clean
# 54581 features across 309734 samples within 1 assay 

HCp.clean <- AddMetaData(HCp.clean, metadata = tmp.meta)

DimPlot(HCp.clean, pt.size = 0.1, group.by = "class_id_label", raster = FALSE, label = TRUE) + NoLegend()

table(HCp.clean@meta.data$class_id_label)
table(HCp.clean@meta.data$subclass_id_label)

DimPlot(HCp.clean, pt.size = 0.1, group.by = "subclass_id_label", raster = FALSE, label = TRUE) + NoLegend()
DimPlot(HCp.clean, pt.size = 0.1, group.by = "supertype_id_label", raster = FALSE, label = TRUE) + NoLegend()


FeaturePlot(HCp.clean, features = c("Snap25"),cols = c("lightgrey","red2"))

FeaturePlot(HCp.clean, features = c("Flt1"),cols = c("lightgrey","red2"))
FeaturePlot(HCp.clean, features = c("Pecam1"),cols = c("lightgrey","red2"))

FeaturePlot(HCp.clean, features = c("Slc17a7"),cols = c("lightgrey","red2"))
FeaturePlot(HCp.clean, features = c("Slc17a8"),cols = c("lightgrey","red2"))
FeaturePlot(HCp.clean, features = c("Gad1"),cols = c("lightgrey","red2"))
FeaturePlot(HCp.clean, features = c("Gad2"),cols = c("lightgrey","red2"))
FeaturePlot(HCp.clean, features = c("Atp1a2"),cols = c("lightgrey","red2"))
FeaturePlot(HCp.clean, features = c("Aqp4"),cols = c("lightgrey","red2"))
FeaturePlot(HCp.clean, features = c("Ctss"),cols = c("lightgrey","red2"))
FeaturePlot(HCp.clean, features = c("Pdgfra"),cols = c("lightgrey","red2"))
FeaturePlot(HCp.clean, features = c("Mbp"),cols = c("lightgrey","red2"))
FeaturePlot(HCp.clean, features = c("Dsp"),cols = c("lightgrey","red2"))
FeaturePlot(HCp.clean, features = c("Hs3st4"),cols = c("lightgrey","red2"))

mcsaveRDS(HCp.clean, file = "./HCp.clean.rna.rds")

#HCp.clean <- mcreadRDS(file = "./HCp.clean.rna.rds")

rm(HCp.clean, HCp, tmp.meta)





############################################################

HYP <- mcreadRDS(file = "./HYP.rna.rds")

tmp.meta <- all.metadata[intersect(rownames(all.metadata), rownames(HYP@meta.data)),c(1:30)]

HYP
# 54588 features across 317448 samples within 1 assay 

HYP.clean <- subset(HYP, cells = rownames(tmp.meta))

HYP.clean
# 54588 features across 313017 samples within 1 assay 

HYP.clean <- AddMetaData(HYP.clean, metadata = tmp.meta)
DimPlot(HYP.clean, pt.size = 0.1, group.by = "class_id_label", raster = FALSE, label = TRUE) + NoLegend()

table(HYP.clean@meta.data$class_id_label)
table(HYP.clean@meta.data$subclass_id_label)

DimPlot(HYP.clean, pt.size = 0.1, group.by = "subclass_id_label", raster = FALSE, label = TRUE) + NoLegend()
DimPlot(HYP.clean, pt.size = 0.1, group.by = "supertype_id_label", raster = FALSE, label = TRUE) + NoLegend()


FeaturePlot(HYP.clean, features = c("Snap25"),cols = c("lightgrey","red2"))

FeaturePlot(HYP.clean, features = c("Flt1"),cols = c("lightgrey","red2"))
FeaturePlot(HYP.clean, features = c("Pecam1"),cols = c("lightgrey","red2"))

FeaturePlot(HYP.clean, features = c("Slc17a7"),cols = c("lightgrey","red2"))
FeaturePlot(HYP.clean, features = c("Slc17a8"),cols = c("lightgrey","red2"))
FeaturePlot(HYP.clean, features = c("Gad1"),cols = c("lightgrey","red2"))
FeaturePlot(HYP.clean, features = c("Gad2"),cols = c("lightgrey","red2"))
FeaturePlot(HYP.clean, features = c("Atp1a2"),cols = c("lightgrey","red2"))
FeaturePlot(HYP.clean, features = c("Aqp4"),cols = c("lightgrey","red2"))
FeaturePlot(HYP.clean, features = c("Ctss"),cols = c("lightgrey","red2"))
FeaturePlot(HYP.clean, features = c("Pdgfra"),cols = c("lightgrey","red2"))
FeaturePlot(HYP.clean, features = c("Mbp"),cols = c("lightgrey","red2"))
FeaturePlot(HYP.clean, features = c("Dsp"),cols = c("lightgrey","red2"))
FeaturePlot(HYP.clean, features = c("Hs3st4"),cols = c("lightgrey","red2"))

mcsaveRDS(HYP.clean, file = "./HYP.clean.rna.rds")

#HYP.clean <- mcreadRDS(file = "./HYP.clean.rna.rds")

rm(HYP.clean, HYP, tmp.meta)


############################################################

NAC <- mcreadRDS(file = "./NAC.rna.rds")

tmp.meta <- all.metadata[intersect(rownames(all.metadata), rownames(NAC@meta.data)),c(1:30)]

NAC
# 54270 features across 250757 samples within 1 assay

NAC.clean <- subset(NAC, cells = rownames(tmp.meta))

NAC.clean
# 54270 features across 235655 samples within 1 assay

NAC.clean <- AddMetaData(NAC.clean, metadata = tmp.meta)
DimPlot(NAC.clean, pt.size = 0.1, group.by = "class_id_label", raster = FALSE, label = TRUE) + NoLegend()

table(NAC.clean@meta.data$class_id_label)
table(NAC.clean@meta.data$subclass_id_label)

DimPlot(NAC.clean, pt.size = 0.1, group.by = "subclass_id_label", raster = FALSE, label = TRUE) + NoLegend()
DimPlot(NAC.clean, pt.size = 0.1, group.by = "supertype_id_label", raster = FALSE, label = TRUE) + NoLegend()


FeaturePlot(NAC.clean, features = c("Snap25"),cols = c("lightgrey","red2"))

FeaturePlot(NAC.clean, features = c("Flt1"),cols = c("lightgrey","red2"))
FeaturePlot(NAC.clean, features = c("Pecam1"),cols = c("lightgrey","red2"))

FeaturePlot(NAC.clean, features = c("Slc17a7"),cols = c("lightgrey","red2"))
FeaturePlot(NAC.clean, features = c("Slc17a8"),cols = c("lightgrey","red2"))
FeaturePlot(NAC.clean, features = c("Gad1"),cols = c("lightgrey","red2"))
FeaturePlot(NAC.clean, features = c("Gad2"),cols = c("lightgrey","red2"))
FeaturePlot(NAC.clean, features = c("Atp1a2"),cols = c("lightgrey","red2"))
FeaturePlot(NAC.clean, features = c("Aqp4"),cols = c("lightgrey","red2"))
FeaturePlot(NAC.clean, features = c("Ctss"),cols = c("lightgrey","red2"))
FeaturePlot(NAC.clean, features = c("Pdgfra"),cols = c("lightgrey","red2"))
FeaturePlot(NAC.clean, features = c("Mbp"),cols = c("lightgrey","red2"))
FeaturePlot(NAC.clean, features = c("Dsp"),cols = c("lightgrey","red2"))
FeaturePlot(NAC.clean, features = c("Hs3st4"),cols = c("lightgrey","red2"))

mcsaveRDS(NAC.clean, file = "./NAC.clean.rna.rds")

#NAC.clean <- mcreadRDS(file = "./NAC.clean.rna.rds")

rm(NAC.clean, NAC, tmp.meta)




############################################################

PFC <- mcreadRDS(file = "./PFC.rna.rds")

tmp.meta <- all.metadata[intersect(rownames(all.metadata), rownames(PFC@meta.data)),c(1:30)]

PFC
# 54270 features across 250757 samples within 1 assay

PFC.clean <- subset(PFC, cells = rownames(tmp.meta))

PFC.clean
# 54270 features across 235655 samples within 1 assay

PFC.clean <- AddMetaData(PFC.clean, metadata = tmp.meta)
DimPlot(PFC.clean, pt.size = 0.1, group.by = "class_id_label", raster = FALSE, label = TRUE) + NoLegend()

table(PFC.clean@meta.data$class_id_label)
table(PFC.clean@meta.data$subclass_id_label)

DimPlot(PFC.clean, pt.size = 0.1, group.by = "subclass_id_label", raster = FALSE, label = TRUE) + NoLegend()
DimPlot(PFC.clean, pt.size = 0.1, group.by = "supertype_id_label", raster = FALSE, label = TRUE) + NoLegend()


FeaturePlot(PFC.clean, features = c("Snap25"),cols = c("lightgrey","red2"))

FeaturePlot(PFC.clean, features = c("Flt1"),cols = c("lightgrey","red2"))
FeaturePlot(PFC.clean, features = c("Pecam1"),cols = c("lightgrey","red2"))

FeaturePlot(PFC.clean, features = c("Slc17a7"),cols = c("lightgrey","red2"))
FeaturePlot(PFC.clean, features = c("Slc17a8"),cols = c("lightgrey","red2"))
FeaturePlot(PFC.clean, features = c("Gad1"),cols = c("lightgrey","red2"))
FeaturePlot(PFC.clean, features = c("Gad2"),cols = c("lightgrey","red2"))
FeaturePlot(PFC.clean, features = c("Atp1a2"),cols = c("lightgrey","red2"))
FeaturePlot(PFC.clean, features = c("Aqp4"),cols = c("lightgrey","red2"))
FeaturePlot(PFC.clean, features = c("Ctss"),cols = c("lightgrey","red2"))
FeaturePlot(PFC.clean, features = c("Pdgfra"),cols = c("lightgrey","red2"))
FeaturePlot(PFC.clean, features = c("Mbp"),cols = c("lightgrey","red2"))
FeaturePlot(PFC.clean, features = c("Dsp"),cols = c("lightgrey","red2"))
FeaturePlot(PFC.clean, features = c("Hs3st4"),cols = c("lightgrey","red2"))

mcsaveRDS(PFC.clean, file = "./PFC.clean.rna.rds")

#PFC.clean <- mcreadRDS(file = "./PFC.clean.rna.rds")

rm(PFC.clean, PFC, tmp.meta)


############################################################

VTA <- mcreadRDS(file = "./VTA_SnR.rna.rds")

tmp.meta <- all.metadata[intersect(rownames(all.metadata), rownames(VTA@meta.data)),c(1:30)]

VTA
# 54300 features across 340000 samples within 1 assay 

VTA.clean <- subset(VTA, cells = rownames(tmp.meta))

VTA.clean
# 54300 features across 336276 samples within 1 assay

VTA.clean <- AddMetaData(VTA.clean, metadata = tmp.meta)
DimPlot(VTA.clean, pt.size = 0.1, group.by = "class_id_label", raster = FALSE, label = TRUE) + NoLegend()

table(VTA.clean@meta.data$class_id_label)
table(VTA.clean@meta.data$subclass_id_label)

DimPlot(VTA.clean, pt.size = 0.1, group.by = "subclass_id_label", raster = FALSE, label = TRUE) + NoLegend()
DimPlot(VTA.clean, pt.size = 0.1, group.by = "supertype_id_label", raster = FALSE, label = TRUE) + NoLegend()


FeaturePlot(VTA.clean, features = c("Snap25"),cols = c("lightgrey","red2"))

FeaturePlot(VTA.clean, features = c("Flt1"),cols = c("lightgrey","red2"))
FeaturePlot(VTA.clean, features = c("Pecam1"),cols = c("lightgrey","red2"))

FeaturePlot(VTA.clean, features = c("Slc17a7"),cols = c("lightgrey","red2"))
FeaturePlot(VTA.clean, features = c("Slc17a8"),cols = c("lightgrey","red2"))
FeaturePlot(VTA.clean, features = c("Gad1"),cols = c("lightgrey","red2"))
FeaturePlot(VTA.clean, features = c("Gad2"),cols = c("lightgrey","red2"))
FeaturePlot(VTA.clean, features = c("Atp1a2"),cols = c("lightgrey","red2"))
FeaturePlot(VTA.clean, features = c("Aqp4"),cols = c("lightgrey","red2"))
FeaturePlot(VTA.clean, features = c("Ctss"),cols = c("lightgrey","red2"))
FeaturePlot(VTA.clean, features = c("Pdgfra"),cols = c("lightgrey","red2"))
FeaturePlot(VTA.clean, features = c("Mbp"),cols = c("lightgrey","red2"))
FeaturePlot(VTA.clean, features = c("Dsp"),cols = c("lightgrey","red2"))
FeaturePlot(VTA.clean, features = c("Hs3st4"),cols = c("lightgrey","red2"))

mcsaveRDS(VTA.clean, file = "./VTA_SnR.clean.rna.rds")

#VTA.clean <- mcreadRDS(file = "./VTA.clean.rna.rds")

rm(VTA.clean, VTA, tmp.meta)


##################################

AMY.clean <- mcreadRDS(file = "./AMY.clean.rna.rds")

AMY.celltype.markers <- read.csv("./AMY.subclass.markers.20240523.txt", header = TRUE)

top1 <- AMY.celltype.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
genes <- unique(top1$gene)
DotPlot(AMY.clean, features=genes, group.by = "subclass_id_label", cols = c("lightgrey","red2")) + RotatedAxis()


# FeatureScatter(AMY.clean, feature1 = "L1UMAP1", feature2 = "L1UMAP2", 
#                group.by = "class_id_label", 
#                cols = c("01 IT-ET Glut" = "red2"),
#                raster = FALSE, pt.size = 0.1)

AMY.metadata <- AMY.clean@meta.data

L1umap <-AMY.metadata[,c("L1UMAP1","L1UMAP2")]

L1umap <- as.matrix(L1umap)

AMY.backup <- AMY.clean

AMY.clean[["L1UMAP"]] <- CreateDimReducObject(embeddings = L1umap, assay = DefaultAssay(AMY.clean), key = "L1UMAP_")

AMY.clean
Idents(AMY.clean) <-"class_id_label"

DimPlot(AMY.clean, reduction = "L1UMAP", 
        group.by = "class_id_label",
        cells.highlight = WhichCells(AMY.clean, idents = "01 IT-ET Glut"),
        cols.highlight = "red2",
        order = TRUE,
        raster = FALSE, 
        pt.size = 0.1,
        sizes.highlight = 0.1)

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(AMY.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./AMY_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(AMY.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(AMY.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}




# Optionally, display all plots (this will display them one by one in a loop)
for (p in umap_plots) {
  print(p)
}


mcsaveRDS(AMY.clean, file = "./AMY.clean.rna.rds")

rm(AMY.clean, AMY.celltype.markers, AMY.metadata, L1umap, umap_plots, top1, cell_groups, output_dir)

rm(VTA, VTA.clean, p, p1, top1)

rm(cell_group, cell_groups, cells_to_highlight,genes, output_dir)

rm(all.metadata, tmp.meta)
############################################

############################################

CPU.clean <- mcreadRDS(file = "./CPU.clean.rna.rds")

CPU.celltype.markers <- read.csv("./CPU.subclass.markers.20240523.txt", header = TRUE)

top1 <- CPU.celltype.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
genes <- unique(top1$gene)
DotPlot(CPU.clean, features=genes, group.by = "subclass_id_label", cols = c("lightgrey","red2")) + RotatedAxis()


# FeatureScatter(CPU.clean, feature1 = "L1UMAP1", feature2 = "L1UMAP2", 
#                group.by = "class_id_label", 
#                cols = c("01 IT-ET Glut" = "red2"),
#                raster = FALSE, pt.size = 0.1)

CPU.metadata <- CPU.clean@meta.data

L1umap <-CPU.metadata[,c("L1UMAP1","L1UMAP2")]

L1umap <- as.matrix(L1umap)

#CPU.backup <- CPU.clean

CPU.clean[["L1UMAP"]] <- CreateDimReducObject(embeddings = L1umap, assay = DefaultAssay(CPU.clean), key = "L1UMAP_")

CPU.clean
Idents(CPU.clean) <-"class_id_label"

# DimPlot(CPU.clean, reduction = "L1UMAP", 
#         group.by = "class_id_label",
#         cells.highlight = WhichCells(CPU.clean, idents = "01 IT-ET Glut"),
#         cols.highlight = "red2",
#         order = TRUE,
#         raster = FALSE, 
#         pt.size = 0.1,
#         sizes.highlight = 0.1)

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(CPU.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./CPU_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(CPU.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(CPU.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}


# # Optionally, display all plots (this will display them one by one in a loop)
# for (p in umap_plots) {
#   print(p)
# }
# 

mcsaveRDS(CPU.clean, file = "./CPU.clean.rna.rds")

rm(CPU.clean, CPU.celltype.markers, CPU.metadata, L1umap, umap_plots, p, top1, cell_group,cell_groups, cells_to_highlight, genes,output_dir)

############################################



############################################

ERC.clean <- mcreadRDS(file = "./ERC.clean.rna.rds")

ERC.celltype.markers <- read.csv("./ERC.subclass.markers.20240523.txt", header = TRUE)

top1 <- ERC.celltype.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
genes <- unique(top1$gene)
DotPlot(ERC.clean, features=genes, group.by = "subclass_id_label", cols = c("lightgrey","red2")) + RotatedAxis()


# FeatureScatter(ERC.clean, feature1 = "L1UMAP1", feature2 = "L1UMAP2", 
#                group.by = "class_id_label", 
#                cols = c("01 IT-ET Glut" = "red2"),
#                raster = FALSE, pt.size = 0.1)

ERC.metadata <- ERC.clean@meta.data

L1umap <-ERC.metadata[,c("L1UMAP1","L1UMAP2")]

L1umap <- as.matrix(L1umap)

#ERC.backup <- ERC.clean

ERC.clean[["L1UMAP"]] <- CreateDimReducObject(embeddings = L1umap, assay = DefaultAssay(ERC.clean), key = "L1UMAP_")

ERC.clean
Idents(ERC.clean) <-"class_id_label"

# DimPlot(ERC.clean, reduction = "L1UMAP", 
#         group.by = "class_id_label",
#         cells.highlight = WhichCells(ERC.clean, idents = "01 IT-ET Glut"),
#         cols.highlight = "red2",
#         order = TRUE,
#         raster = FALSE, 
#         pt.size = 0.1,
#         sizes.highlight = 0.1)

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(ERC.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./ERC_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(ERC.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(ERC.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}


# # Optionally, display all plots (this will display them one by one in a loop)
# for (p in umap_plots) {
#   print(p)
# }
# 

mcsaveRDS(ERC.clean, file = "./ERC.clean.rna.rds")

rm(ERC.clean, ERC.celltype.markers, ERC.metadata, L1umap, umap_plots, p, top1, cell_group,cell_groups, cells_to_highlight, genes,output_dir)

############################################



############################################

HCa.clean <- mcreadRDS(file = "./HCa.clean.rna.rds")

HCa.celltype.markers <- read.csv("./HIP.subclass.markers.20240523.txt", header = TRUE)

top1 <- HCa.celltype.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
genes <- unique(top1$gene)

markerplot <- DotPlot(HCa.clean, features=genes, group.by = "subclass_id_label", cols = c("lightgrey","red2")) + RotatedAxis()

pdf("./20240609.HCa.PT.subclass.marker.pdf", width=16, height=10)
print(markerplot)
dev.off()




# FeatureScatter(HCa.clean, feature1 = "L1UMAP1", feature2 = "L1UMAP2", 
#                group.by = "class_id_label", 
#                cols = c("01 IT-ET Glut" = "red2"),
#                raster = FALSE, pt.size = 0.1)

HCa.metadata <- HCa.clean@meta.data

L1umap <-HCa.metadata[,c("L1UMAP1","L1UMAP2")]

L1umap <- as.matrix(L1umap)

#HCa.backup <- HCa.clean

HCa.clean[["L1UMAP"]] <- CreateDimReducObject(embeddings = L1umap, assay = DefaultAssay(HCa.clean), key = "L1UMAP_")

HCa.clean
Idents(HCa.clean) <-"class_id_label"

# DimPlot(HCa.clean, reduction = "L1UMAP", 
#         group.by = "class_id_label",
#         cells.highlight = WhichCells(HCa.clean, idents = "01 IT-ET Glut"),
#         cols.highlight = "red2",
#         order = TRUE,
#         raster = FALSE, 
#         pt.size = 0.1,
#         sizes.highlight = 0.1)

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(HCa.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./HCa_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(HCa.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(HCa.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}


# # Optionally, display all plots (this will display them one by one in a loop)
# for (p in umap_plots) {
#   print(p)
# }
# 

mcsaveRDS(HCa.clean, file = "./HCa.clean.rna.rds")

rm(HCa.clean, HCa.celltype.markers, HCa.metadata, L1umap, umap_plots, p, top1, cell_group,cell_groups, cells_to_highlight, genes,output_dir)

############################################



############################################

HCp.clean <- mcreadRDS(file = "./HCp.clean.rna.rds")

HCp.celltype.markers <- read.csv("./HIP.subclass.markers.20240523.txt", header = TRUE)

top1 <- HCp.celltype.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
genes <- unique(top1$gene)

markerplot <- DotPlot(HCp.clean, features=genes, group.by = "subclass_id_label", cols = c("lightgrey","red2")) + RotatedAxis()

pdf("./20240609.HCp.PT.subclass.marker.pdf", width=16, height=10)
print(markerplot)
dev.off()



# FeatureScatter(HCp.clean, feature1 = "L1UMAP1", feature2 = "L1UMAP2", 
#                group.by = "class_id_label", 
#                cols = c("01 IT-ET Glut" = "red2"),
#                raster = FALSE, pt.size = 0.1)

HCp.metadata <- HCp.clean@meta.data

L1umap <-HCp.metadata[,c("L1UMAP1","L1UMAP2")]

L1umap <- as.matrix(L1umap)

#HCp.backup <- HCp.clean

HCp.clean[["L1UMAP"]] <- CreateDimReducObject(embeddings = L1umap, assay = DefaultAssay(HCp.clean), key = "L1UMAP_")

HCp.clean
Idents(HCp.clean) <-"class_id_label"

# DimPlot(HCp.clean, reduction = "L1UMAP", 
#         group.by = "class_id_label",
#         cells.highlight = WhichCells(HCp.clean, idents = "01 IT-ET Glut"),
#         cols.highlight = "red2",
#         order = TRUE,
#         raster = FALSE, 
#         pt.size = 0.1,
#         sizes.highlight = 0.1)

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(HCp.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./HCp_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(HCp.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(HCp.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}


# # Optionally, display all plots (this will display them one by one in a loop)
# for (p in umap_plots) {
#   print(p)
# }
# 

mcsaveRDS(HCp.clean, file = "./HCp.clean.rna.rds")

rm(HCp.clean, HCp.celltype.markers, HCp.metadata, L1umap, umap_plots, p, top1, cell_group,cell_groups, cells_to_highlight, genes,output_dir)

############################################



############################################

HYP.clean <- mcreadRDS(file = "./HYP.clean.rna.rds")

HYP.celltype.markers <- read.csv("./HYP.subclass.markers.20240523.txt", header = TRUE)

top1 <- HYP.celltype.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
genes <- unique(top1$gene)

markerplot <- DotPlot(HYP.clean, features=genes, group.by = "subclass_id_label", cols = c("lightgrey","red2")) + RotatedAxis()

pdf("./20240609.HYP.PT.subclass.marker.pdf", width=16, height=10)
print(markerplot)
dev.off()



# FeatureScatter(HYP.clean, feature1 = "L1UMAP1", feature2 = "L1UMAP2", 
#                group.by = "class_id_label", 
#                cols = c("01 IT-ET Glut" = "red2"),
#                raster = FALSE, pt.size = 0.1)

HYP.metadata <- HYP.clean@meta.data

L1umap <-HYP.metadata[,c("L1UMAP1","L1UMAP2")]

L1umap <- as.matrix(L1umap)

#HYP.backup <- HYP.clean

HYP.clean[["L1UMAP"]] <- CreateDimReducObject(embeddings = L1umap, assay = DefaultAssay(HYP.clean), key = "L1UMAP_")

HYP.clean
Idents(HYP.clean) <-"class_id_label"

# DimPlot(HYP.clean, reduction = "L1UMAP", 
#         group.by = "class_id_label",
#         cells.highlight = WhichCells(HYP.clean, idents = "01 IT-ET Glut"),
#         cols.highlight = "red2",
#         order = TRUE,
#         raster = FALSE, 
#         pt.size = 0.1,
#         sizes.highlight = 0.1)

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(HYP.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./HYP_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(HYP.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(HYP.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}


# # Optionally, display all plots (this will display them one by one in a loop)
# for (p in umap_plots) {
#   print(p)
# }
# 

mcsaveRDS(HYP.clean, file = "./HYP.clean.rna.rds")

rm(HYP.clean, HYP.celltype.markers, HYP.metadata, L1umap, umap_plots, p, top1, cell_group,cell_groups, cells_to_highlight, genes,output_dir)

############################################





############################################

NAC.clean <- mcreadRDS(file = "./NAC.clean.rna.rds")

NAC.celltype.markers <- read.csv("./NAC.subclass.markers.20240523.txt", header = TRUE)

top1 <- NAC.celltype.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
genes <- unique(top1$gene)


markerplot <- DotPlot(NAC.clean, features=genes, group.by = "subclass_id_label", cols = c("lightgrey","red2")) + RotatedAxis()

pdf("./20240609.NAC.PT.subclass.marker.pdf", width=16, height=10)
print(markerplot)
dev.off()



# FeatureScatter(NAC.clean, feature1 = "L1UMAP1", feature2 = "L1UMAP2", 
#                group.by = "class_id_label", 
#                cols = c("01 IT-ET Glut" = "red2"),
#                raster = FALSE, pt.size = 0.1)

NAC.metadata <- NAC.clean@meta.data

L1umap <-NAC.metadata[,c("L1UMAP1","L1UMAP2")]

L1umap <- as.matrix(L1umap)

#NAC.backup <- NAC.clean

NAC.clean[["L1UMAP"]] <- CreateDimReducObject(embeddings = L1umap, assay = DefaultAssay(NAC.clean), key = "L1UMAP_")

NAC.clean
Idents(NAC.clean) <-"class_id_label"

# DimPlot(NAC.clean, reduction = "L1UMAP", 
#         group.by = "class_id_label",
#         cells.highlight = WhichCells(NAC.clean, idents = "01 IT-ET Glut"),
#         cols.highlight = "red2",
#         order = TRUE,
#         raster = FALSE, 
#         pt.size = 0.1,
#         sizes.highlight = 0.1)

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(NAC.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./NAC_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(NAC.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(NAC.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}


# # Optionally, display all plots (this will display them one by one in a loop)
# for (p in umap_plots) {
#   print(p)
# }
# 

mcsaveRDS(NAC.clean, file = "./NAC.clean.rna.rds")

rm(NAC.clean, NAC.celltype.markers, NAC.metadata, L1umap, umap_plots, p, top1, cell_group,cell_groups, cells_to_highlight, genes,output_dir)

############################################






############################################

PFC.clean <- mcreadRDS(file = "./PFC.clean.rna.rds")

PFC.celltype.markers <- read.csv("./PFC.subclass.markers.20240523.txt", header = TRUE)

top1 <- PFC.celltype.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
genes <- unique(top1$gene)


markerplot <- DotPlot(PFC.clean, features=genes, group.by = "subclass_id_label", cols = c("lightgrey","red2")) + RotatedAxis()

pdf("./20240609.PFC.PT.subclass.marker.pdf", width=16, height=10)
print(markerplot)
dev.off()


# FeatureScatter(PFC.clean, feature1 = "L1UMAP1", feature2 = "L1UMAP2", 
#                group.by = "class_id_label", 
#                cols = c("01 IT-ET Glut" = "red2"),
#                raster = FALSE, pt.size = 0.1)

PFC.metadata <- PFC.clean@meta.data

L1umap <-PFC.metadata[,c("L1UMAP1","L1UMAP2")]

L1umap <- as.matrix(L1umap)

#PFC.backup <- PFC.clean

PFC.clean[["L1UMAP"]] <- CreateDimReducObject(embeddings = L1umap, assay = DefaultAssay(PFC.clean), key = "L1UMAP_")

PFC.clean
Idents(PFC.clean) <-"class_id_label"

# DimPlot(PFC.clean, reduction = "L1UMAP", 
#         group.by = "class_id_label",
#         cells.highlight = WhichCells(PFC.clean, idents = "01 IT-ET Glut"),
#         cols.highlight = "red2",
#         order = TRUE,
#         raster = FALSE, 
#         pt.size = 0.1,
#         sizes.highlight = 0.1)

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(PFC.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./PFC_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(PFC.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(PFC.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}


# # Optionally, display all plots (this will display them one by one in a loop)
# for (p in umap_plots) {
#   print(p)
# }
# 

mcsaveRDS(PFC.clean, file = "./PFC.clean.rna.rds")

rm(PFC.clean, PFC.celltype.markers, PFC.metadata, L1umap, umap_plots, p, top1, cell_group,cell_groups, cells_to_highlight, genes,output_dir)

############################################




############################################

VTA_SnR.clean <- mcreadRDS(file = "./VTA_SnR.clean.rna.rds")

VTA_SnR.celltype.markers <- read.csv("./VTA.subclass.markers.20240523.txt", header = TRUE)

top1 <- VTA_SnR.celltype.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
genes <- unique(top1$gene)

markerplot <- DotPlot(VTA_SnR.clean, features=genes, group.by = "subclass_id_label", cols = c("lightgrey","red2")) + RotatedAxis()

pdf("./20240609.VTA_SnR.PT.subclass.marker.pdf", width=16, height=10)
print(markerplot)
dev.off()



# FeatureScatter(VTA_SnR.clean, feature1 = "L1UMAP1", feature2 = "L1UMAP2", 
#                group.by = "class_id_label", 
#                cols = c("01 IT-ET Glut" = "red2"),
#                raster = FALSE, pt.size = 0.1)

VTA_SnR.metadata <- VTA_SnR.clean@meta.data

L1umap <-VTA_SnR.metadata[,c("L1UMAP1","L1UMAP2")]

L1umap <- as.matrix(L1umap)

#VTA_SnR.backup <- VTA_SnR.clean

VTA_SnR.clean[["L1UMAP"]] <- CreateDimReducObject(embeddings = L1umap, assay = DefaultAssay(VTA_SnR.clean), key = "L1UMAP_")

VTA_SnR.clean
Idents(VTA_SnR.clean) <-"class_id_label"

# DimPlot(VTA_SnR.clean, reduction = "L1UMAP", 
#         group.by = "class_id_label",
#         cells.highlight = WhichCells(VTA_SnR.clean, idents = "01 IT-ET Glut"),
#         cols.highlight = "red2",
#         order = TRUE,
#         raster = FALSE, 
#         pt.size = 0.1,
#         sizes.highlight = 0.1)

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(VTA_SnR.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./VTA_SnR_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(VTA_SnR.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}


# # Optionally, display all plots (this will display them one by one in a loop)
# for (p in umap_plots) {
#   print(p)
# }
# 

mcsaveRDS(VTA_SnR.clean, file = "./VTA_SnR.clean.rna.rds")

rm(VTA_SnR.clean, VTA_SnR.celltype.markers, VTA_SnR.metadata, L1umap, umap_plots, p, top1, cell_group,cell_groups, cells_to_highlight, genes,output_dir)

############################################


AMY.clean <- mcreadRDS(file = "./AMY.clean.rna.rds")

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(AMY.clean@meta.data$class_id_label)

# > cell_groups
# [1] "01 IT-ET Glut"   "09 CNU-LGE GABA" "11 CNU-HYa GABA" "13 CNU-HYa Glut" "30 Astro-Epen"   "31 OPC-Oligo"    "33 Vascular"    
# [8] "08 CNU-MGE GABA" "34 Immune"       "04 DG-IMN Glut"  "07 CTX-MGE GABA" "06 CTX-CGE GABA" "03 OB-CR Glut"   "05 OB-IMN GABA" 



# Initialize an empty list to store the plots
umap_plots <- list()

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(AMY.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(AMY.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
        ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p

}

AMY.clean[["umap_filter"]] <- "remove"

table(AMY.clean@meta.data$umap_filter)

#################
umap.plot1 <- umap_plots$`01 IT-ET Glut`
Idents(AMY.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(AMY.clean, idents = "01 IT-ET Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(AMY.clean) <- "umap_filter"
Idents(AMY.clean, cells = select.cells) <- "keep"
AMY.clean[["umap_filter"]] <- Idents(AMY.clean)

table(AMY.clean@meta.data$umap_filter)
###########################

check <- subset(x = AMY.clean, subset = class_id_label =="01 IT-ET Glut")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", split.by = "subclass_id_label", ncol = 4)

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(AMY.clean@meta.data$subclass_id_label)
DimPlot(AMY.clean, group.by = "subclass_id_label", reduction = "L1UMAP")


#################
umap.plot1 <- umap_plots$`09 CNU-LGE GABA`
Idents(AMY.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(AMY.clean, idents = "09 CNU-LGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(AMY.clean) <- "umap_filter"
Idents(AMY.clean, cells = select.cells) <- "keep"
AMY.clean[["umap_filter"]] <- Idents(AMY.clean)
###########################

check <- subset(x = AMY.clean, subset = class_id_label =="09 CNU-LGE GABA")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", split.by = "subclass_id_label", ncol = 3)

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(AMY.clean@meta.data$subclass_id_label)
DimPlot(AMY.clean, group.by = "subclass_id_label", reduction = "L1UMAP")



#################
umap.plot1 <- umap_plots$`11 CNU-HYa GABA`
Idents(AMY.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(AMY.clean, idents = "11 CNU-HYa GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(AMY.clean) <- "umap_filter"
Idents(AMY.clean, cells = select.cells) <- "keep"
AMY.clean[["umap_filter"]] <- Idents(AMY.clean)

###########################

check <- subset(x = AMY.clean, subset = class_id_label =="11 CNU-HYa GABA")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(AMY.clean@meta.data$subclass_id_label)
DimPlot(AMY.clean, group.by = "subclass_id_label", reduction = "L1UMAP")


#################
umap.plot1 <- umap_plots$`13 CNU-HYa Glut`
Idents(AMY.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(AMY.clean, idents = "13 CNU-HYa Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(AMY.clean) <- "umap_filter"
Idents(AMY.clean, cells = select.cells) <- "keep"
AMY.clean[["umap_filter"]] <- Idents(AMY.clean)


###########################

check <- subset(x = AMY.clean, subset = class_id_label =="13 CNU-HYa Glut")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(AMY.clean@meta.data$subclass_id_label)
DimPlot(AMY.clean, group.by = "subclass_id_label", reduction = "L1UMAP")






#################
umap.plot1 <- umap_plots$`30 Astro-Epen`
Idents(AMY.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(AMY.clean, idents = "30 Astro-Epen")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(AMY.clean) <- "umap_filter"
Idents(AMY.clean, cells = select.cells) <- "keep"
AMY.clean[["umap_filter"]] <- Idents(AMY.clean)


#################
umap.plot1 <- umap_plots$`31 OPC-Oligo`
Idents(AMY.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(AMY.clean, idents = "31 OPC-Oligo")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(AMY.clean) <- "umap_filter"
Idents(AMY.clean, cells = select.cells) <- "keep"
AMY.clean[["umap_filter"]] <- Idents(AMY.clean)

#################
umap.plot1 <- umap_plots$`33 Vascular`
Idents(AMY.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(AMY.clean, idents = "33 Vascular")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(AMY.clean) <- "umap_filter"
Idents(AMY.clean, cells = select.cells) <- "keep"
AMY.clean[["umap_filter"]] <- Idents(AMY.clean)

#################
umap.plot1 <- umap_plots$`08 CNU-MGE GABA`
Idents(AMY.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(AMY.clean, idents = "08 CNU-MGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(AMY.clean) <- "umap_filter"
Idents(AMY.clean, cells = select.cells) <- "keep"
AMY.clean[["umap_filter"]] <- Idents(AMY.clean)

#################
umap.plot1 <- umap_plots$`34 Immune`
Idents(AMY.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(AMY.clean, idents = "34 Immune")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(AMY.clean) <- "umap_filter"
Idents(AMY.clean, cells = select.cells) <- "keep"
AMY.clean[["umap_filter"]] <- Idents(AMY.clean)


###########################

check <- subset(x = AMY.clean, subset = class_id_label =="34 Immune")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(AMY.clean@meta.data$subclass_id_label)
DimPlot(AMY.clean, group.by = "subclass_id_label", reduction = "L1UMAP")


#################
umap.plot1 <- umap_plots$`04 DG-IMN Glut`
Idents(AMY.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(AMY.clean, idents = "04 DG-IMN Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(AMY.clean) <- "umap_filter"
Idents(AMY.clean, cells = select.cells) <- "keep"
AMY.clean[["umap_filter"]] <- Idents(AMY.clean)


###########################

check <- subset(x = AMY.clean, subset = class_id_label =="04 DG-IMN Glut")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")
DimPlot(check, group.by = "supertype_id_label", reduction = "L1UMAP")

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(AMY.clean@meta.data$subclass_id_label)
DimPlot(AMY.clean, group.by = "subclass_id_label", reduction = "L1UMAP")




#################
umap.plot1 <- umap_plots$`07 CTX-MGE GABA`
Idents(AMY.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(AMY.clean, idents = "07 CTX-MGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(AMY.clean) <- "umap_filter"
Idents(AMY.clean, cells = select.cells) <- "keep"
AMY.clean[["umap_filter"]] <- Idents(AMY.clean)


###########################

check <- subset(x = AMY.clean, subset = class_id_label =="07 CTX-MGE GABA")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")
DimPlot(check, group.by = "supertype_id_label", reduction = "L1UMAP")

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(AMY.clean@meta.data$subclass_id_label)
DimPlot(AMY.clean, group.by = "subclass_id_label", reduction = "L1UMAP")




#################
umap.plot1 <- umap_plots$`06 CTX-CGE GABA`
Idents(AMY.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(AMY.clean, idents = "06 CTX-CGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(AMY.clean) <- "umap_filter"
Idents(AMY.clean, cells = select.cells) <- "keep"
AMY.clean[["umap_filter"]] <- Idents(AMY.clean)

#################
umap.plot1 <- umap_plots$`03 OB-CR Glut`
Idents(AMY.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(AMY.clean, idents = "03 OB-CR Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(AMY.clean) <- "umap_filter"
Idents(AMY.clean, cells = select.cells) <- "keep"
AMY.clean[["umap_filter"]] <- Idents(AMY.clean)


###########################

check <- subset(x = AMY.clean, subset = class_id_label =="03 OB-CR Glut")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")
DimPlot(check, group.by = "supertype_id_label", reduction = "L1UMAP")

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(AMY.clean@meta.data$subclass_id_label)
DimPlot(AMY.clean, group.by = "subclass_id_label", reduction = "L1UMAP")



#################
umap.plot1 <- umap_plots$`05 OB-IMN GABA`
Idents(AMY.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(AMY.clean, idents = "05 OB-IMN GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(AMY.clean) <- "umap_filter"
Idents(AMY.clean, cells = select.cells) <- "keep"
AMY.clean[["umap_filter"]] <- Idents(AMY.clean)

###########################

check <- subset(x = AMY.clean, subset = class_id_label =="05 OB-IMN GABA")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")
DimPlot(check, group.by = "supertype_id_label", reduction = "L1UMAP")

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(AMY.clean@meta.data$subclass_id_label)
DimPlot(AMY.clean, group.by = "subclass_id_label", reduction = "L1UMAP")

###################

#################
umap.plot1 <- umap_plots$`02 NP-CT-L6b Glut`
Idents(AMY.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(AMY.clean, idents = "02 NP-CT-L6b Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(AMY.clean) <- "umap_filter"
Idents(AMY.clean, cells = select.cells) <- "keep"
AMY.clean[["umap_filter"]] <- Idents(AMY.clean)

###########################

check <- subset(x = AMY.clean, subset = class_id_label =="02 NP-CT-L6b Glut")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")
DimPlot(check, group.by = "supertype_id_label", reduction = "L1UMAP")

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(AMY.clean@meta.data$subclass_id_label)
DimPlot(AMY.clean, group.by = "subclass_id_label", reduction = "L1UMAP")

###################




######################
DimPlot(AMY.clean, group.by = "umap_filter", reduction = "L1UMAP")

table(AMY.clean@meta.data$umap_filter, AMY.clean@meta.data$class_id_label)

table(AMY.clean@meta.data$umap_filter)
# keep    remove 
# 243123  15680 
######################


############################################
mcsaveRDS(AMY.clean, file = "./AMY.clean.rna.rds")

Idents(AMY.clean) <-"umap_filter"

AMY.clean <- subset(AMY.clean, idents = c("keep"))
Idents(AMY.clean) <-"class_id_label"

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(AMY.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./AMY_clean_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(AMY.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(AMY.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}



rm(AMY.clean, umap_plots, p, cell_group,cell_groups, cells_to_highlight, output_dir)

############################################

CPU.clean <- mcreadRDS(file = "./CPU.clean.rna.rds")

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(CPU.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(CPU.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(CPU.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
}

CPU.clean[["umap_filter"]] <- "remove"

table(CPU.clean@meta.data$umap_filter)

#################
umap.plot1 <- umap_plots$`09 CNU-LGE GABA`
Idents(CPU.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(CPU.clean, idents = "09 CNU-LGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(CPU.clean) <- "umap_filter"
Idents(CPU.clean, cells = select.cells) <- "keep"
CPU.clean[["umap_filter"]] <- Idents(CPU.clean)

#################
umap.plot1 <- umap_plots$`30 Astro-Epen`
Idents(CPU.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(CPU.clean, idents = "30 Astro-Epen")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(CPU.clean) <- "umap_filter"
Idents(CPU.clean, cells = select.cells) <- "keep"
CPU.clean[["umap_filter"]] <- Idents(CPU.clean)

#################
umap.plot1 <- umap_plots$`31 OPC-Oligo`
Idents(CPU.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(CPU.clean, idents = "31 OPC-Oligo")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(CPU.clean) <- "umap_filter"
Idents(CPU.clean, cells = select.cells) <- "keep"
CPU.clean[["umap_filter"]] <- Idents(CPU.clean)

#################
umap.plot1 <- umap_plots$`02 NP-CT-L6b Glut`
Idents(CPU.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(CPU.clean, idents = "02 NP-CT-L6b Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(CPU.clean) <- "umap_filter"
Idents(CPU.clean, cells = select.cells) <- "keep"
CPU.clean[["umap_filter"]] <- Idents(CPU.clean)


#################
umap.plot1 <- umap_plots$`34 Immune`
Idents(CPU.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(CPU.clean, idents = "34 Immune")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(CPU.clean) <- "umap_filter"
Idents(CPU.clean, cells = select.cells) <- "keep"
CPU.clean[["umap_filter"]] <- Idents(CPU.clean)


#################
umap.plot1 <- umap_plots$`08 CNU-MGE GABA`
Idents(CPU.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(CPU.clean, idents = "08 CNU-MGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(CPU.clean) <- "umap_filter"
Idents(CPU.clean, cells = select.cells) <- "keep"
CPU.clean[["umap_filter"]] <- Idents(CPU.clean)


#################
umap.plot1 <- umap_plots$`05 OB-IMN GABA`
Idents(CPU.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(CPU.clean, idents = "05 OB-IMN GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(CPU.clean) <- "umap_filter"
Idents(CPU.clean, cells = select.cells) <- "keep"
CPU.clean[["umap_filter"]] <- Idents(CPU.clean)

#################
umap.plot1 <- umap_plots$`06 CTX-CGE GABA`
Idents(CPU.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(CPU.clean, idents = "06 CTX-CGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(CPU.clean) <- "umap_filter"
Idents(CPU.clean, cells = select.cells) <- "keep"
CPU.clean[["umap_filter"]] <- Idents(CPU.clean)

#################
umap.plot1 <- umap_plots$`33 Vascular`
Idents(CPU.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(CPU.clean, idents = "33 Vascular")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(CPU.clean) <- "umap_filter"
Idents(CPU.clean, cells = select.cells) <- "keep"
CPU.clean[["umap_filter"]] <- Idents(CPU.clean)


#################
umap.plot1 <- umap_plots$`01 IT-ET Glut`
Idents(CPU.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(CPU.clean, idents = "01 IT-ET Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(CPU.clean) <- "umap_filter"
Idents(CPU.clean, cells = select.cells) <- "keep"
CPU.clean[["umap_filter"]] <- Idents(CPU.clean)


#################
umap.plot1 <- umap_plots$`07 CTX-MGE GABA`
Idents(CPU.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(CPU.clean, idents = "07 CTX-MGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(CPU.clean) <- "umap_filter"
Idents(CPU.clean, cells = select.cells) <- "keep"
CPU.clean[["umap_filter"]] <- Idents(CPU.clean)


#################
umap.plot1 <- umap_plots$`04 DG-IMN Glut`
# Very sparse, delete entire cluster

#################
umap.plot1 <- umap_plots$`11 CNU-HYa GABA`
Idents(CPU.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(CPU.clean, idents = "11 CNU-HYa GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(CPU.clean) <- "umap_filter"
Idents(CPU.clean, cells = select.cells) <- "keep"
CPU.clean[["umap_filter"]] <- Idents(CPU.clean)


###########################

check <- subset(x = CPU.clean, subset = class_id_label =="05 OB-IMN GABA")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(CPU.clean@meta.data$subclass_id_label)
DimPlot(CPU.clean, group.by = "subclass_id_label", reduction = "L1UMAP")

######################

table(CPU.clean@meta.data$umap_filter)
# keep   remove 
# 337171  12295 
######################



############################################
mcsaveRDS(CPU.clean, file = "./CPU.clean.rna.rds")

Idents(CPU.clean) <-"umap_filter"

CPU.clean <- subset(CPU.clean, idents = c("keep"))
Idents(CPU.clean) <-"class_id_label"

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(CPU.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./CPU_clean_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(CPU.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(CPU.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}



rm(CPU.clean, CPU.celltype.markers, CPU.metadata, L1umap, umap_plots, p, top1, cell_group,cell_groups, cells_to_highlight, genes,output_dir)

############################################



############################################

ERC.clean <- mcreadRDS(file = "./ERC.clean.rna.rds")

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(ERC.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(ERC.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(ERC.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
}

ERC.clean[["umap_filter"]] <- "remove"

table(ERC.clean@meta.data$umap_filter)

#################
umap.plot1 <- umap_plots$`01 IT-ET Glut`
Idents(ERC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(ERC.clean, idents = "01 IT-ET Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(ERC.clean) <- "umap_filter"
Idents(ERC.clean, cells = select.cells) <- "keep"
ERC.clean[["umap_filter"]] <- Idents(ERC.clean)


#################
umap.plot1 <- umap_plots$`33 Vascular`
Idents(ERC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(ERC.clean, idents = "33 Vascular")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(ERC.clean) <- "umap_filter"
Idents(ERC.clean, cells = select.cells) <- "keep"
ERC.clean[["umap_filter"]] <- Idents(ERC.clean)



#################
umap.plot1 <- umap_plots$`30 Astro-Epen`
Idents(ERC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(ERC.clean, idents = "30 Astro-Epen")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(ERC.clean) <- "umap_filter"
Idents(ERC.clean, cells = select.cells) <- "keep"
ERC.clean[["umap_filter"]] <- Idents(ERC.clean)


#################
umap.plot1 <- umap_plots$`02 NP-CT-L6b Glut`
Idents(ERC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(ERC.clean, idents = "02 NP-CT-L6b Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(ERC.clean) <- "umap_filter"
Idents(ERC.clean, cells = select.cells) <- "keep"
ERC.clean[["umap_filter"]] <- Idents(ERC.clean)


#################
umap.plot1 <- umap_plots$`31 OPC-Oligo`
Idents(ERC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(ERC.clean, idents = "31 OPC-Oligo")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(ERC.clean) <- "umap_filter"
Idents(ERC.clean, cells = select.cells) <- "keep"
ERC.clean[["umap_filter"]] <- Idents(ERC.clean)

#################
umap.plot1 <- umap_plots$`07 CTX-MGE GABA`
Idents(ERC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(ERC.clean, idents = "07 CTX-MGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(ERC.clean) <- "umap_filter"
Idents(ERC.clean, cells = select.cells) <- "keep"
ERC.clean[["umap_filter"]] <- Idents(ERC.clean)

#################
umap.plot1 <- umap_plots$`06 CTX-CGE GABA`
Idents(ERC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(ERC.clean, idents = "06 CTX-CGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(ERC.clean) <- "umap_filter"
Idents(ERC.clean, cells = select.cells) <- "keep"
ERC.clean[["umap_filter"]] <- Idents(ERC.clean)


#################
umap.plot1 <- umap_plots$`08 CNU-MGE GABA`
Idents(ERC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(ERC.clean, idents = "08 CNU-MGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(ERC.clean) <- "umap_filter"
Idents(ERC.clean, cells = select.cells) <- "keep"
ERC.clean[["umap_filter"]] <- Idents(ERC.clean)


#################
umap.plot1 <- umap_plots$`34 Immune`
Idents(ERC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(ERC.clean, idents = "34 Immune")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(ERC.clean) <- "umap_filter"
Idents(ERC.clean, cells = select.cells) <- "keep"
ERC.clean[["umap_filter"]] <- Idents(ERC.clean)

#################
umap.plot1 <- umap_plots$`04 DG-IMN Glut`
Idents(ERC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(ERC.clean, idents = "04 DG-IMN Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(ERC.clean) <- "umap_filter"
Idents(ERC.clean, cells = select.cells) <- "keep"
ERC.clean[["umap_filter"]] <- Idents(ERC.clean)

#################
umap.plot1 <- umap_plots$`05 OB-IMN GABA`
Idents(ERC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(ERC.clean, idents = "05 OB-IMN GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(ERC.clean) <- "umap_filter"
Idents(ERC.clean, cells = select.cells) <- "keep"
ERC.clean[["umap_filter"]] <- Idents(ERC.clean)


#################
umap.plot1 <- umap_plots$`03 OB-CR Glut`
Idents(ERC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(ERC.clean, idents = "03 OB-CR Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(ERC.clean) <- "umap_filter"
Idents(ERC.clean, cells = select.cells) <- "keep"
ERC.clean[["umap_filter"]] <- Idents(ERC.clean)





###########################

check <- subset(x = ERC.clean, subset = class_id_label =="04 DG-IMN Glut")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(ERC.clean@meta.data$subclass_id_label)
DimPlot(ERC.clean, group.by = "subclass_id_label", reduction = "L1UMAP")

######################

table(ERC.clean@meta.data$umap_filter)
# keep   remove 
# 307986  21363 
######################



############################################
mcsaveRDS(ERC.clean, file = "./ERC.clean.rna.rds")

Idents(ERC.clean) <-"umap_filter"

ERC.clean <- subset(ERC.clean, idents = c("keep"))
Idents(ERC.clean) <-"class_id_label"

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(ERC.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./ERC_clean_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(ERC.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(ERC.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}



rm(ERC.clean, ERC.celltype.markers, ERC.metadata, L1umap, umap_plots, p, top1, cell_group,cell_groups, cells_to_highlight, genes,output_dir)

############################################


############################################

HCa.clean <- mcreadRDS(file = "./HCa.clean.rna.rds")

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(HCa.clean@meta.data$class_id_label)


# Initialize an empty list to store the plots
umap_plots <- list()
Idents(HCa.clean) <- "class_id_label"
# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(HCa.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(HCa.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
}

HCa.clean[["umap_filter"]] <- "remove"

table(HCa.clean@meta.data$umap_filter)

#################
umap.plot1 <- umap_plots$`04 DG-IMN Glut`
Idents(HCa.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCa.clean, idents = "04 DG-IMN Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCa.clean) <- "umap_filter"
Idents(HCa.clean, cells = select.cells) <- "keep"
HCa.clean[["umap_filter"]] <- Idents(HCa.clean)


#################
umap.plot1 <- umap_plots$`34 Immune`
Idents(HCa.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCa.clean, idents = "34 Immune")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCa.clean) <- "umap_filter"
Idents(HCa.clean, cells = select.cells) <- "keep"
HCa.clean[["umap_filter"]] <- Idents(HCa.clean)

#################
umap.plot1 <- umap_plots$`30 Astro-Epen`
Idents(HCa.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCa.clean, idents = "30 Astro-Epen")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCa.clean) <- "umap_filter"
Idents(HCa.clean, cells = select.cells) <- "keep"
HCa.clean[["umap_filter"]] <- Idents(HCa.clean)


#################
umap.plot1 <- umap_plots$`01 IT-ET Glut`
Idents(HCa.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCa.clean, idents = "01 IT-ET Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCa.clean) <- "umap_filter"
Idents(HCa.clean, cells = select.cells) <- "keep"
HCa.clean[["umap_filter"]] <- Idents(HCa.clean)

#################
umap.plot1 <- umap_plots$`31 OPC-Oligo`
Idents(HCa.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCa.clean, idents = "31 OPC-Oligo")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCa.clean) <- "umap_filter"
Idents(HCa.clean, cells = select.cells) <- "keep"
HCa.clean[["umap_filter"]] <- Idents(HCa.clean)


#################
umap.plot1 <- umap_plots$`07 CTX-MGE GABA`
Idents(HCa.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCa.clean, idents = "07 CTX-MGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCa.clean) <- "umap_filter"
Idents(HCa.clean, cells = select.cells) <- "keep"
HCa.clean[["umap_filter"]] <- Idents(HCa.clean)


#################
umap.plot1 <- umap_plots$`03 OB-CR Glut`
Idents(HCa.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCa.clean, idents = "03 OB-CR Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCa.clean) <- "umap_filter"
Idents(HCa.clean, cells = select.cells) <- "keep"
HCa.clean[["umap_filter"]] <- Idents(HCa.clean)



######################
mcsaveRDS(HCa.clean, file = "./HCa.clean.rna.rds")

rm(HCa.clean)

AMY.clean <- mcreadRDS(file = "./AMY.clean.rna.rds")

Idents(AMY.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(AMY.clean, idents = "03 OB-CR Glut")
p <- DimPlot(AMY.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle("03 OB-CR Glut")
  

#################
umap.plot1 <- p
Idents(AMY.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(AMY.clean, idents = "03 OB-CR Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(AMY.clean) <- "umap_filter"
Idents(AMY.clean, cells = select.cells) <- "remove"
AMY.clean[["umap_filter"]] <- Idents(AMY.clean)

table(AMY.clean@meta.data$umap_filter)


Idents(AMY.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(AMY.clean, idents = "03 OB-CR Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(AMY.clean) <- "umap_filter"
Idents(AMY.clean, cells = select.cells) <- "keep"
AMY.clean[["umap_filter"]] <- Idents(AMY.clean)


############################################
mcsaveRDS(AMY.clean, file = "./AMY.clean.rna.rds")

Idents(AMY.clean) <-"umap_filter"

AMY.clean <- subset(AMY.clean, idents = c("keep"))
Idents(AMY.clean) <-"class_id_label"

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(AMY.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./AMY_clean_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(AMY.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(AMY.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}



rm(AMY.clean, umap_plots, p, cell_group,cell_groups, cells_to_highlight, output_dir)

# I am pausing here! - 20240610 12:57AM



#################
umap.plot1 <- umap_plots$`03 OB-CR Glut`
Idents(HCa.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCa.clean, idents = "03 OB-CR Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCa.clean) <- "umap_filter"
Idents(HCa.clean, cells = select.cells) <- "keep"
HCa.clean[["umap_filter"]] <- Idents(HCa.clean)


#################
umap.plot1 <- umap_plots$`06 CTX-CGE GABA`
Idents(HCa.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCa.clean, idents = "06 CTX-CGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCa.clean) <- "umap_filter"
Idents(HCa.clean, cells = select.cells) <- "keep"
HCa.clean[["umap_filter"]] <- Idents(HCa.clean)

#################
umap.plot1 <- umap_plots$`33 Vascular`
Idents(HCa.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCa.clean, idents = "33 Vascular")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCa.clean) <- "umap_filter"
Idents(HCa.clean, cells = select.cells) <- "keep"
HCa.clean[["umap_filter"]] <- Idents(HCa.clean)

#################
umap.plot1 <- umap_plots$`02 NP-CT-L6b Glut`
Idents(HCa.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCa.clean, idents = "02 NP-CT-L6b Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCa.clean) <- "umap_filter"
Idents(HCa.clean, cells = select.cells) <- "keep"
HCa.clean[["umap_filter"]] <- Idents(HCa.clean)


#################
umap.plot1 <- umap_plots$`05 OB-IMN GABA`
Idents(HCa.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCa.clean, idents = "05 OB-IMN GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCa.clean) <- "umap_filter"
Idents(HCa.clean, cells = select.cells) <- "keep"
HCa.clean[["umap_filter"]] <- Idents(HCa.clean)


###########################

check <- subset(x = HCa.clean, subset = class_id_label =="03 OB-CR Glut")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(HCa.clean@meta.data$subclass_id_label)
DimPlot(HCa.clean, group.by = "subclass_id_label", reduction = "L1UMAP")
DimPlot(HCa.clean, group.by = "umap_filter", reduction = "L1UMAP")
DimPlot(HCa.clean, group.by = "class_id_label", reduction = "L1UMAP")





######################

table(HCa.clean@meta.data$umap_filter)
# keep    remove 
# 402045   9842 
######################



############################################
mcsaveRDS(HCa.clean, file = "./HCa.clean.rna.rds")

Idents(HCa.clean) <-"umap_filter"

HCa.clean <- subset(HCa.clean, idents = c("keep"))
Idents(HCa.clean) <-"class_id_label"

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(HCa.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./HCa_clean_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(HCa.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(HCa.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}



rm(HCa.clean, HCa.celltype.markers, HCa.metadata, L1umap, umap_plots, p, top1, cell_group,cell_groups, cells_to_highlight, genes,output_dir)

############################################

HCp.clean <- mcreadRDS(file = "./HCp.clean.rna.rds")

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(HCp.clean@meta.data$class_id_label)


# Initialize an empty list to store the plots
umap_plots <- list()

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(HCp.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(HCp.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
}

HCp.clean[["umap_filter"]] <- "remove"

table(HCp.clean@meta.data$umap_filter)

#################
umap.plot1 <- umap_plots$`31 OPC-Oligo`
Idents(HCp.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCp.clean, idents = "31 OPC-Oligo")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCp.clean) <- "umap_filter"
Idents(HCp.clean, cells = select.cells) <- "keep"
HCp.clean[["umap_filter"]] <- Idents(HCp.clean)


#################
umap.plot1 <- umap_plots$`04 DG-IMN Glut`
Idents(HCp.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCp.clean, idents = "04 DG-IMN Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCp.clean) <- "umap_filter"
Idents(HCp.clean, cells = select.cells) <- "keep"
HCp.clean[["umap_filter"]] <- Idents(HCp.clean)

#################
umap.plot1 <- umap_plots$`01 IT-ET Glut`
Idents(HCp.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCp.clean, idents = "01 IT-ET Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCp.clean) <- "umap_filter"
Idents(HCp.clean, cells = select.cells) <- "keep"
HCp.clean[["umap_filter"]] <- Idents(HCp.clean)

#################
umap.plot1 <- umap_plots$`03 OB-CR Glut`
Idents(HCp.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCp.clean, idents = "03 OB-CR Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCp.clean) <- "umap_filter"
Idents(HCp.clean, cells = select.cells) <- "keep"
HCp.clean[["umap_filter"]] <- Idents(HCp.clean)

#################
umap.plot1 <- umap_plots$`34 Immune`
Idents(HCp.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCp.clean, idents = "34 Immune")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCp.clean) <- "umap_filter"
Idents(HCp.clean, cells = select.cells) <- "keep"
HCp.clean[["umap_filter"]] <- Idents(HCp.clean)


#################
umap.plot1 <- umap_plots$`30 Astro-Epen`
Idents(HCp.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCp.clean, idents = "30 Astro-Epen")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCp.clean) <- "umap_filter"
Idents(HCp.clean, cells = select.cells) <- "keep"
HCp.clean[["umap_filter"]] <- Idents(HCp.clean)


#################
umap.plot1 <- umap_plots$`07 CTX-MGE GABA`
Idents(HCp.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCp.clean, idents = "07 CTX-MGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCp.clean) <- "umap_filter"
Idents(HCp.clean, cells = select.cells) <- "keep"
HCp.clean[["umap_filter"]] <- Idents(HCp.clean)


#################
umap.plot1 <- umap_plots$`06 CTX-CGE GABA`
Idents(HCp.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCp.clean, idents = "06 CTX-CGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCp.clean) <- "umap_filter"
Idents(HCp.clean, cells = select.cells) <- "keep"
HCp.clean[["umap_filter"]] <- Idents(HCp.clean)



#################
umap.plot1 <- umap_plots$`02 NP-CT-L6b Glut`
Idents(HCp.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCp.clean, idents = "02 NP-CT-L6b Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCp.clean) <- "umap_filter"
Idents(HCp.clean, cells = select.cells) <- "keep"
HCp.clean[["umap_filter"]] <- Idents(HCp.clean)


#################
umap.plot1 <- umap_plots$`33 Vascular`
Idents(HCp.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCp.clean, idents = "33 Vascular")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCp.clean) <- "umap_filter"
Idents(HCp.clean, cells = select.cells) <- "keep"
HCp.clean[["umap_filter"]] <- Idents(HCp.clean)



#################
umap.plot1 <- umap_plots$`05 OB-IMN GABA`
Idents(HCp.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HCp.clean, idents = "05 OB-IMN GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HCp.clean) <- "umap_filter"
Idents(HCp.clean, cells = select.cells) <- "keep"
HCp.clean[["umap_filter"]] <- Idents(HCp.clean)



###########################

check <- subset(x = HCp.clean, subset = class_id_label =="33 Vascular")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(HCp.clean@meta.data$subclass_id_label)
DimPlot(HCp.clean, group.by = "subclass_id_label", reduction = "L1UMAP")

######################

table(HCp.clean@meta.data$umap_filter)
# keep    remove 
# 298025  11709 
######################



############################################
mcsaveRDS(HCp.clean, file = "./HCp.clean.rna.rds")

Idents(HCp.clean) <-"umap_filter"

HCp.clean <- subset(HCp.clean, idents = c("keep"))
Idents(HCp.clean) <-"class_id_label"

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(HCp.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./HCp_clean_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(HCp.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(HCp.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}



rm(HCp.clean, umap_plots, p, cell_group,cell_groups, cells_to_highlight,output_dir)
gc()
############################################


############################################

HYP.clean <- mcreadRDS(file = "./HYP.clean.rna.rds")

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(HYP.clean@meta.data$class_id_label)


# Initialize an empty list to store the plots
umap_plots <- list()

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(HYP.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(HYP.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
}

HYP.clean[["umap_filter"]] <- "remove"

table(HYP.clean@meta.data$umap_filter)

#################
umap.plot1 <- umap_plots$`30 Astro-Epen`
Idents(HYP.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HYP.clean, idents = "30 Astro-Epen")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HYP.clean) <- "umap_filter"
Idents(HYP.clean, cells = select.cells) <- "keep"
HYP.clean[["umap_filter"]] <- Idents(HYP.clean)


#################
umap.plot1 <- umap_plots$`15 HY Gnrh1 Glut`
Idents(HYP.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HYP.clean, idents = "15 HY Gnrh1 Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HYP.clean) <- "umap_filter"
Idents(HYP.clean, cells = select.cells) <- "keep"
HYP.clean[["umap_filter"]] <- Idents(HYP.clean)

#################
umap.plot1 <- umap_plots$`12 HY GABA`
Idents(HYP.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HYP.clean, idents = "12 HY GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HYP.clean) <- "umap_filter"
Idents(HYP.clean, cells = select.cells) <- "keep"
HYP.clean[["umap_filter"]] <- Idents(HYP.clean)


#################
umap.plot1 <- umap_plots$`11 CNU-HYa GABA`
Idents(HYP.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HYP.clean, idents = "11 CNU-HYa GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HYP.clean) <- "umap_filter"
Idents(HYP.clean, cells = select.cells) <- "keep"
HYP.clean[["umap_filter"]] <- Idents(HYP.clean)


#################
umap.plot1 <- umap_plots$`31 OPC-Oligo`
Idents(HYP.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HYP.clean, idents = "31 OPC-Oligo")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HYP.clean) <- "umap_filter"
Idents(HYP.clean, cells = select.cells) <- "keep"
HYP.clean[["umap_filter"]] <- Idents(HYP.clean)

#################
umap.plot1 <- umap_plots$`14 HY Glut`
Idents(HYP.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HYP.clean, idents = "14 HY Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HYP.clean) <- "umap_filter"
Idents(HYP.clean, cells = select.cells) <- "keep"
HYP.clean[["umap_filter"]] <- Idents(HYP.clean)

#################
umap.plot1 <- umap_plots$`18 TH Glut`
Idents(HYP.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HYP.clean, idents = "18 TH Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HYP.clean) <- "umap_filter"
Idents(HYP.clean, cells = select.cells) <- "keep"
HYP.clean[["umap_filter"]] <- Idents(HYP.clean)


#################
umap.plot1 <- umap_plots$`13 CNU-HYa Glut`
Idents(HYP.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HYP.clean, idents = "13 CNU-HYa Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HYP.clean) <- "umap_filter"
Idents(HYP.clean, cells = select.cells) <- "keep"
HYP.clean[["umap_filter"]] <- Idents(HYP.clean)

#################
umap.plot1 <- umap_plots$`34 Immune`
Idents(HYP.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HYP.clean, idents = "34 Immune")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HYP.clean) <- "umap_filter"
Idents(HYP.clean, cells = select.cells) <- "keep"
HYP.clean[["umap_filter"]] <- Idents(HYP.clean)


#################
umap.plot1 <- umap_plots$`33 Vascular`
Idents(HYP.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HYP.clean, idents = "33 Vascular")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HYP.clean) <- "umap_filter"
Idents(HYP.clean, cells = select.cells) <- "keep"
HYP.clean[["umap_filter"]] <- Idents(HYP.clean)

#################
umap.plot1 <- umap_plots$`08 CNU-MGE GABA`
Idents(HYP.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HYP.clean, idents = "08 CNU-MGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HYP.clean) <- "umap_filter"
Idents(HYP.clean, cells = select.cells) <- "keep"
HYP.clean[["umap_filter"]] <- Idents(HYP.clean)


#################
umap.plot1 <- umap_plots$`19 MB Glut`
Idents(HYP.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HYP.clean, idents = "19 MB Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HYP.clean) <- "umap_filter"
Idents(HYP.clean, cells = select.cells) <- "keep"
HYP.clean[["umap_filter"]] <- Idents(HYP.clean)



#################
umap.plot1 <- umap_plots$`04 DG-IMN Glut`
Idents(HYP.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HYP.clean, idents = "04 DG-IMN Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HYP.clean) <- "umap_filter"
Idents(HYP.clean, cells = select.cells) <- "keep"
HYP.clean[["umap_filter"]] <- Idents(HYP.clean)


#################
umap.plot1 <- umap_plots$`20 MB GABA`
Idents(HYP.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HYP.clean, idents = "20 MB GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HYP.clean) <- "umap_filter"
Idents(HYP.clean, cells = select.cells) <- "keep"
HYP.clean[["umap_filter"]] <- Idents(HYP.clean)



#################
umap.plot1 <- umap_plots$`09 CNU-LGE GABA`
Idents(HYP.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HYP.clean, idents = "09 CNU-LGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HYP.clean) <- "umap_filter"
Idents(HYP.clean, cells = select.cells) <- "keep"
HYP.clean[["umap_filter"]] <- Idents(HYP.clean)

#################
umap.plot1 <- umap_plots$`05 OB-IMN GABA`
Idents(HYP.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HYP.clean, idents = "05 OB-IMN GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HYP.clean) <- "umap_filter"
Idents(HYP.clean, cells = select.cells) <- "keep"
HYP.clean[["umap_filter"]] <- Idents(HYP.clean)


#################
umap.plot1 <- umap_plots$`16 HY MM Glut`
Idents(HYP.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(HYP.clean, idents = "16 HY MM Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(HYP.clean) <- "umap_filter"
Idents(HYP.clean, cells = select.cells) <- "keep"
HYP.clean[["umap_filter"]] <- Idents(HYP.clean)




###########################

check <- subset(x = HYP.clean, subset = class_id_label =="05 OB-IMN GABA")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(HYP.clean@meta.data$subclass_id_label)
DimPlot(HYP.clean, group.by = "subclass_id_label", reduction = "L1UMAP")

######################

table(HYP.clean@meta.data$umap_filter)
# keep    remove 
# 308988   4029 
######################



############################################
mcsaveRDS(HYP.clean, file = "./HYP.clean.rna.rds")

Idents(HYP.clean) <-"umap_filter"

HYP.clean <- subset(HYP.clean, idents = c("keep"))
Idents(HYP.clean) <-"class_id_label"

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(HYP.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./HYP_clean_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(HYP.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(HYP.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}



rm(HYP.clean, umap_plots, p, cell_group,cell_groups, cells_to_highlight,output_dir)
gc()

############################################



############################################

NAC.clean <- mcreadRDS(file = "./NAC.clean.rna.rds")

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(NAC.clean@meta.data$class_id_label)


# Initialize an empty list to store the plots
umap_plots <- list()

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(NAC.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(NAC.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
}

NAC.clean[["umap_filter"]] <- "remove"

table(NAC.clean@meta.data$umap_filter)

#################
umap.plot1 <- umap_plots$`09 CNU-LGE GABA`
Idents(NAC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(NAC.clean, idents = "09 CNU-LGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(NAC.clean) <- "umap_filter"
Idents(NAC.clean, cells = select.cells) <- "keep"
NAC.clean[["umap_filter"]] <- Idents(NAC.clean)

table(NAC.clean@meta.data$subclass_id_label, NAC.clean@meta.data$umap_filter)


#################
umap.plot1 <- umap_plots$`31 OPC-Oligo`
Idents(NAC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(NAC.clean, idents = "31 OPC-Oligo")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(NAC.clean) <- "umap_filter"
Idents(NAC.clean, cells = select.cells) <- "keep"
NAC.clean[["umap_filter"]] <- Idents(NAC.clean)

#################
umap.plot1 <- umap_plots$`30 Astro-Epen`
Idents(NAC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(NAC.clean, idents = "30 Astro-Epen")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(NAC.clean) <- "umap_filter"
Idents(NAC.clean, cells = select.cells) <- "keep"
NAC.clean[["umap_filter"]] <- Idents(NAC.clean)


#################
umap.plot1 <- umap_plots$`05 OB-IMN GABA`
Idents(NAC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(NAC.clean, idents = "05 OB-IMN GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(NAC.clean) <- "umap_filter"
Idents(NAC.clean, cells = select.cells) <- "keep"
NAC.clean[["umap_filter"]] <- Idents(NAC.clean)

#################
umap.plot1 <- umap_plots$`34 Immune`
Idents(NAC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(NAC.clean, idents = "34 Immune")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(NAC.clean) <- "umap_filter"
Idents(NAC.clean, cells = select.cells) <- "keep"
NAC.clean[["umap_filter"]] <- Idents(NAC.clean)


#################
umap.plot1 <- umap_plots$`01 IT-ET Glut`
Idents(NAC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(NAC.clean, idents = "01 IT-ET Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(NAC.clean) <- "umap_filter"
Idents(NAC.clean, cells = select.cells) <- "keep"
NAC.clean[["umap_filter"]] <- Idents(NAC.clean)


#################
umap.plot1 <- umap_plots$`11 CNU-HYa GABA`
Idents(NAC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(NAC.clean, idents = "11 CNU-HYa GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(NAC.clean) <- "umap_filter"
Idents(NAC.clean, cells = select.cells) <- "keep"
NAC.clean[["umap_filter"]] <- Idents(NAC.clean)


#################
umap.plot1 <- umap_plots$`08 CNU-MGE GABA`
Idents(NAC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(NAC.clean, idents = "08 CNU-MGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(NAC.clean) <- "umap_filter"
Idents(NAC.clean, cells = select.cells) <- "keep"
NAC.clean[["umap_filter"]] <- Idents(NAC.clean)



#################
umap.plot1 <- umap_plots$`33 Vascular`
Idents(NAC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(NAC.clean, idents = "33 Vascular")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(NAC.clean) <- "umap_filter"
Idents(NAC.clean, cells = select.cells) <- "keep"
NAC.clean[["umap_filter"]] <- Idents(NAC.clean)



#################
umap.plot1 <- umap_plots$`02 NP-CT-L6b Glut`
Idents(NAC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(NAC.clean, idents = "02 NP-CT-L6b Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(NAC.clean) <- "umap_filter"
Idents(NAC.clean, cells = select.cells) <- "keep"
NAC.clean[["umap_filter"]] <- Idents(NAC.clean)


#################
umap.plot1 <- umap_plots$`04 DG-IMN Glut`
Idents(NAC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(NAC.clean, idents = "04 DG-IMN Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(NAC.clean) <- "umap_filter"
Idents(NAC.clean, cells = select.cells) <- "keep"
NAC.clean[["umap_filter"]] <- Idents(NAC.clean)



#################
umap.plot1 <- umap_plots$`07 CTX-MGE GABA`
Idents(NAC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(NAC.clean, idents = "07 CTX-MGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(NAC.clean) <- "umap_filter"
Idents(NAC.clean, cells = select.cells) <- "keep"
NAC.clean[["umap_filter"]] <- Idents(NAC.clean)

#################
umap.plot1 <- umap_plots$`06 CTX-CGE GABA`
Idents(NAC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(NAC.clean, idents = "06 CTX-CGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(NAC.clean) <- "umap_filter"
Idents(NAC.clean, cells = select.cells) <- "keep"
NAC.clean[["umap_filter"]] <- Idents(NAC.clean)


#################
umap.plot1 <- umap_plots$`03 OB-CR Glut`
umap.plot1
# Sparse, and not in the right umap space. delete entire cluster

#################
umap.plot1 <- umap_plots$`13 CNU-HYa Glut`
Idents(NAC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(NAC.clean, idents = "13 CNU-HYa Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(NAC.clean) <- "umap_filter"
Idents(NAC.clean, cells = select.cells) <- "keep"
NAC.clean[["umap_filter"]] <- Idents(NAC.clean)


#################
umap.plot1 <- umap_plots$`10 LSX GABA`
Idents(NAC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(NAC.clean, idents = "10 LSX GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(NAC.clean) <- "umap_filter"
Idents(NAC.clean, cells = select.cells) <- "keep"
NAC.clean[["umap_filter"]] <- Idents(NAC.clean)



###########################

check <- subset(x = NAC.clean, subset = class_id_label =="03 OB-CR Glut")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(NAC.clean@meta.data$subclass_id_label)
DimPlot(NAC.clean, group.by = "subclass_id_label", reduction = "L1UMAP")
DimPlot(NAC.clean, group.by = "umap_filter", reduction = "L1UMAP")

######################

table(NAC.clean@meta.data$umap_filter)
# keep    remove 
# 223472  12183 
######################



############################################
mcsaveRDS(NAC.clean, file = "./NAC.clean.rna.rds")

Idents(NAC.clean) <-"umap_filter"

NAC.clean <- subset(NAC.clean, idents = c("keep"))
Idents(NAC.clean) <-"class_id_label"

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(NAC.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./NAC_clean_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(NAC.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(NAC.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}



rm(NAC.clean, umap_plots, p, cell_group,cell_groups, cells_to_highlight,output_dir)
gc()

############################################



############################################

PFC.clean <- mcreadRDS(file = "./PFC.clean.rna.rds")

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(PFC.clean@meta.data$class_id_label)


# Initialize an empty list to store the plots
umap_plots <- list()

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(PFC.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(PFC.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
}

PFC.clean[["umap_filter"]] <- "remove"

table(PFC.clean@meta.data$umap_filter)

#################
umap.plot1 <- umap_plots$`01 IT-ET Glut`
Idents(PFC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(PFC.clean, idents = "01 IT-ET Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(PFC.clean) <- "umap_filter"
Idents(PFC.clean, cells = select.cells) <- "keep"
PFC.clean[["umap_filter"]] <- Idents(PFC.clean)

#################
umap.plot1 <- umap_plots$`02 NP-CT-L6b Glut`
Idents(PFC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(PFC.clean, idents = "02 NP-CT-L6b Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(PFC.clean) <- "umap_filter"
Idents(PFC.clean, cells = select.cells) <- "keep"
PFC.clean[["umap_filter"]] <- Idents(PFC.clean)

#################
umap.plot1 <- umap_plots$`33 Vascular`
Idents(PFC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(PFC.clean, idents = "33 Vascular")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(PFC.clean) <- "umap_filter"
Idents(PFC.clean, cells = select.cells) <- "keep"
PFC.clean[["umap_filter"]] <- Idents(PFC.clean)



#################
umap.plot1 <- umap_plots$`30 Astro-Epen`
Idents(PFC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(PFC.clean, idents = "30 Astro-Epen")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(PFC.clean) <- "umap_filter"
Idents(PFC.clean, cells = select.cells) <- "keep"
PFC.clean[["umap_filter"]] <- Idents(PFC.clean)



#################
umap.plot1 <- umap_plots$`07 CTX-MGE GABA`
Idents(PFC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(PFC.clean, idents = "07 CTX-MGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(PFC.clean) <- "umap_filter"
Idents(PFC.clean, cells = select.cells) <- "keep"
PFC.clean[["umap_filter"]] <- Idents(PFC.clean)

#################
umap.plot1 <- umap_plots$`34 Immune`
Idents(PFC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(PFC.clean, idents = "34 Immune")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(PFC.clean) <- "umap_filter"
Idents(PFC.clean, cells = select.cells) <- "keep"
PFC.clean[["umap_filter"]] <- Idents(PFC.clean)


#################
umap.plot1 <- umap_plots$`31 OPC-Oligo`
Idents(PFC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(PFC.clean, idents = "31 OPC-Oligo")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(PFC.clean) <- "umap_filter"
Idents(PFC.clean, cells = select.cells) <- "keep"
PFC.clean[["umap_filter"]] <- Idents(PFC.clean)


#################
umap.plot1 <- umap_plots$`06 CTX-CGE GABA`
Idents(PFC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(PFC.clean, idents = "06 CTX-CGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(PFC.clean) <- "umap_filter"
Idents(PFC.clean, cells = select.cells) <- "keep"
PFC.clean[["umap_filter"]] <- Idents(PFC.clean)


#################
umap.plot1 <- umap_plots$`04 DG-IMN Glut`
Idents(PFC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(PFC.clean, idents = "04 DG-IMN Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(PFC.clean) <- "umap_filter"
Idents(PFC.clean, cells = select.cells) <- "keep"
PFC.clean[["umap_filter"]] <- Idents(PFC.clean)


#################
umap.plot1 <- umap_plots$`05 OB-IMN GABA`
Idents(PFC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(PFC.clean, idents = "05 OB-IMN GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(PFC.clean) <- "umap_filter"
Idents(PFC.clean, cells = select.cells) <- "keep"
PFC.clean[["umap_filter"]] <- Idents(PFC.clean)

#################
umap.plot1 <- umap_plots$`08 CNU-MGE GABA`
Idents(PFC.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(PFC.clean, idents = "08 CNU-MGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(PFC.clean) <- "umap_filter"
Idents(PFC.clean, cells = select.cells) <- "keep"
PFC.clean[["umap_filter"]] <- Idents(PFC.clean)


###########################

check <- subset(x = PFC.clean, subset = class_id_label =="08 CNU-MGE GABA")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(PFC.clean@meta.data$class_id_label, PFC.clean@meta.data$umap_filter)
DimPlot(PFC.clean, group.by = "subclass_id_label", reduction = "L1UMAP")
DimPlot(PFC.clean, group.by = "umap_filter", reduction = "L1UMAP")

######################

table(PFC.clean@meta.data$umap_filter)
# keep  remove 
# 148336  11610 
######################



############################################
mcsaveRDS(PFC.clean, file = "./PFC.clean.rna.rds")

Idents(PFC.clean) <-"umap_filter"

PFC.clean <- subset(PFC.clean, idents = c("keep"))
Idents(PFC.clean) <-"class_id_label"

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(PFC.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./PFC_clean_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(PFC.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(PFC.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}



rm(PFC.clean, umap_plots, p, cell_group,cell_groups, cells_to_highlight,output_dir)
gc()

############################################

############################################

VTA_SnR.clean <- mcreadRDS(file = "./VTA_SnR.clean.rna.rds")

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(VTA_SnR.clean@meta.data$class_id_label)


# Initialize an empty list to store the plots
umap_plots <- list()

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(VTA_SnR.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
}

VTA_SnR.clean[["umap_filter"]] <- "remove"

table(VTA_SnR.clean@meta.data$umap_filter)

#################
umap.plot1 <- umap_plots$`31 OPC-Oligo`
Idents(VTA_SnR.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = "31 OPC-Oligo")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(VTA_SnR.clean) <- "umap_filter"
Idents(VTA_SnR.clean, cells = select.cells) <- "keep"
VTA_SnR.clean[["umap_filter"]] <- Idents(VTA_SnR.clean)

#################
umap.plot1 <- umap_plots$`20 MB GABA`
Idents(VTA_SnR.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = "20 MB GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(VTA_SnR.clean) <- "umap_filter"
Idents(VTA_SnR.clean, cells = select.cells) <- "keep"
VTA_SnR.clean[["umap_filter"]] <- Idents(VTA_SnR.clean)

#################
umap.plot1 <- umap_plots$`30 Astro-Epen`
Idents(VTA_SnR.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = "30 Astro-Epen")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(VTA_SnR.clean) <- "umap_filter"
Idents(VTA_SnR.clean, cells = select.cells) <- "keep"
VTA_SnR.clean[["umap_filter"]] <- Idents(VTA_SnR.clean)

#################
umap.plot1 <- umap_plots$`34 Immune`
Idents(VTA_SnR.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = "34 Immune")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(VTA_SnR.clean) <- "umap_filter"
Idents(VTA_SnR.clean, cells = select.cells) <- "keep"
VTA_SnR.clean[["umap_filter"]] <- Idents(VTA_SnR.clean)



#################
umap.plot1 <- umap_plots$`21 MB Dopa`
Idents(VTA_SnR.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = "21 MB Dopa")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(VTA_SnR.clean) <- "umap_filter"
Idents(VTA_SnR.clean, cells = select.cells) <- "keep"
VTA_SnR.clean[["umap_filter"]] <- Idents(VTA_SnR.clean)


#################
umap.plot1 <- umap_plots$`19 MB Glut`
Idents(VTA_SnR.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = "19 MB Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(VTA_SnR.clean) <- "umap_filter"
Idents(VTA_SnR.clean, cells = select.cells) <- "keep"
VTA_SnR.clean[["umap_filter"]] <- Idents(VTA_SnR.clean)

#################
umap.plot1 <- umap_plots$`33 Vascular`
Idents(VTA_SnR.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = "33 Vascular")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(VTA_SnR.clean) <- "umap_filter"
Idents(VTA_SnR.clean, cells = select.cells) <- "keep"
VTA_SnR.clean[["umap_filter"]] <- Idents(VTA_SnR.clean)


#################
umap.plot1 <- umap_plots$`14 HY Glut`
Idents(VTA_SnR.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = "14 HY Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(VTA_SnR.clean) <- "umap_filter"
Idents(VTA_SnR.clean, cells = select.cells) <- "keep"
VTA_SnR.clean[["umap_filter"]] <- Idents(VTA_SnR.clean)


#################
umap.plot1 <- umap_plots$`26 P GABA`
Idents(VTA_SnR.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = "26 P GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(VTA_SnR.clean) <- "umap_filter"
Idents(VTA_SnR.clean, cells = select.cells) <- "keep"
VTA_SnR.clean[["umap_filter"]] <- Idents(VTA_SnR.clean)


################
umap.plot1 <- umap_plots$`22 MB-HB Sero`
Idents(VTA_SnR.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = "22 MB-HB Sero")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(VTA_SnR.clean) <- "umap_filter"
Idents(VTA_SnR.clean, cells = select.cells) <- "keep"
VTA_SnR.clean[["umap_filter"]] <- Idents(VTA_SnR.clean)



################
umap.plot1 <- umap_plots$`04 DG-IMN Glut`
Idents(VTA_SnR.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = "04 DG-IMN Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(VTA_SnR.clean) <- "umap_filter"
Idents(VTA_SnR.clean, cells = select.cells) <- "keep"
VTA_SnR.clean[["umap_filter"]] <- Idents(VTA_SnR.clean)


################
umap.plot1 <- umap_plots$`01 IT-ET Glut`
Idents(VTA_SnR.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = "01 IT-ET Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(VTA_SnR.clean) <- "umap_filter"
Idents(VTA_SnR.clean, cells = select.cells) <- "keep"
VTA_SnR.clean[["umap_filter"]] <- Idents(VTA_SnR.clean)


################
umap.plot1 <- umap_plots$`12 HY GABA`
Idents(VTA_SnR.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = "12 HY GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(VTA_SnR.clean) <- "umap_filter"
Idents(VTA_SnR.clean, cells = select.cells) <- "keep"
VTA_SnR.clean[["umap_filter"]] <- Idents(VTA_SnR.clean)


################
umap.plot1 <- umap_plots$`23 P Glut`
Idents(VTA_SnR.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = "23 P Glut")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(VTA_SnR.clean) <- "umap_filter"
Idents(VTA_SnR.clean, cells = select.cells) <- "keep"
VTA_SnR.clean[["umap_filter"]] <- Idents(VTA_SnR.clean)



################
umap.plot1 <- umap_plots$`05 OB-IMN GABA`
Idents(VTA_SnR.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = "05 OB-IMN GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(VTA_SnR.clean) <- "umap_filter"
Idents(VTA_SnR.clean, cells = select.cells) <- "keep"
VTA_SnR.clean[["umap_filter"]] <- Idents(VTA_SnR.clean)


################
umap.plot1 <- umap_plots$`07 CTX-MGE GABA`
Idents(VTA_SnR.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = "07 CTX-MGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(VTA_SnR.clean) <- "umap_filter"
Idents(VTA_SnR.clean, cells = select.cells) <- "keep"
VTA_SnR.clean[["umap_filter"]] <- Idents(VTA_SnR.clean)

################
umap.plot1 <- umap_plots$`06 CTX-CGE GABA`
Idents(VTA_SnR.clean) <- "class_id_label"
cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = "06 CTX-CGE GABA")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(VTA_SnR.clean) <- "umap_filter"
Idents(VTA_SnR.clean, cells = select.cells) <- "keep"
VTA_SnR.clean[["umap_filter"]] <- Idents(VTA_SnR.clean)



###########################

check <- subset(x = VTA_SnR.clean, subset = class_id_label =="06 CTX-CGE GABA")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")
DimPlot(check, group.by = "supertype_id_label", reduction = "L1UMAP")

FeaturePlot(check, features = c("Th", "Slc6a3"), reduction = "L1UMAP", pt.size = 0.1)

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(VTA_SnR.clean@meta.data$subclass_id_label, VTA_SnR.clean@meta.data$umap_filter)
DimPlot(VTA_SnR.clean, group.by = "subclass_id_label", reduction = "L1UMAP")

DimPlot(VTA_SnR.clean, group.by = "umap_filter", reduction = "L1UMAP")

######################

table(VTA_SnR.clean@meta.data$umap_filter)
# keep   remove 
# 319124  17152 
######################



############################################
mcsaveRDS(VTA_SnR.clean, file = "./VTA_SnR.clean.rna.rds")

Idents(VTA_SnR.clean) <-"umap_filter"

VTA_SnR.clean <- subset(VTA_SnR.clean, idents = c("keep"))
Idents(VTA_SnR.clean) <-"class_id_label"

###########################################
# A loop function to generate highlighted umap plots

# Get unique cell groups from the "class_id_label" column
cell_groups <- unique(VTA_SnR.clean@meta.data$class_id_label)

# Initialize an empty list to store the plots
umap_plots <- list()

# Set output directory for PDF files
output_dir <- "./VTA_SnR_clean_ambPT_class_id_label/"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each cell group and generate a UMAP plot
for (cell_group in cell_groups) {
  # Highlight cells in the current group
  cells_to_highlight <- WhichCells(VTA_SnR.clean, idents = cell_group)
  
  # Generate the UMAP plot
  p <- DimPlot(VTA_SnR.clean, reduction = "L1UMAP", 
               group.by = "class_id_label",
               cells.highlight = cells_to_highlight,
               cols.highlight = "red2",
               order = TRUE,
               raster = FALSE, 
               pt.size = 0.1,
               sizes.highlight = 0.1) +
    ggtitle(cell_group)
  
  # Store the plot in the list
  umap_plots[[cell_group]] <- p
  
  # Save the plot to a PDF file
  pdf(file = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  
}


# Alternatively, save each plot to a file
for (cell_group in names(umap_plots)) {
  ggsave(filename = paste0(output_dir, "/", cell_group, "_UMAP_clean_plot.png"),
         plot = umap_plots[[cell_group]],
         width = 16, height = 12, dpi = 400)
}



rm(VTA_SnR.clean, umap_plots, p, cell_group,cell_groups, cells_to_highlight,output_dir)
gc()

############################################
HCa.clean <- mcreadRDS(file = "./HCa.clean.rna.rds")
DimPlot(HCa.clean, group.by = "class_id_label", reduction = "L1UMAP")
DimPlot(HCa.clean, group.by = "umap_filter", reduction = "L1UMAP")
table(HCa.clean@meta.data$umap_filter)
############################################


AMY.clean <- mcreadRDS(file = "./AMY.clean.rna.rds")
AMY.metadata <- AMY.clean@meta.data
rm(AMY.clean);gc()

CPU.clean <- mcreadRDS(file = "./CPU.clean.rna.rds")
CPU.metadata <- CPU.clean@meta.data
rm(CPU.clean);gc()

ERC.clean <- mcreadRDS(file = "./ERC.clean.rna.rds")
ERC.metadata <- ERC.clean@meta.data
rm(ERC.clean);gc()

HCa.clean <- mcreadRDS(file = "./HCa.clean.rna.rds")
HCa.metadata <- HCa.clean@meta.data
rm(HCa.clean);gc()

HCp.clean <- mcreadRDS(file = "./HCp.clean.rna.rds")
HCp.metadata <- HCp.clean@meta.data
rm(HCp.clean);gc()

HYP.clean <- mcreadRDS(file = "./HYP.clean.rna.rds")
HYP.metadata <- HYP.clean@meta.data
rm(HYP.clean);gc()

NAC.clean <- mcreadRDS(file = "./NAC.clean.rna.rds")
NAC.metadata <- NAC.clean@meta.data
rm(NAC.clean);gc()

PFC.clean <- mcreadRDS(file = "./PFC.clean.rna.rds")
PFC.metadata <- PFC.clean@meta.data
rm(PFC.clean);gc()

VTA_SnR.clean <- mcreadRDS(file = "./VTA_SnR.clean.rna.rds")
VTA_SnR.metadata <- VTA_SnR.clean@meta.data
rm(VTA_SnR.clean);gc()


allcellmetadata <- dplyr::bind_rows(AMY.metadata, CPU.metadata, ERC.metadata, 
                                    HCa.metadata, HCp.metadata, HYP.metadata,
                                    NAC.metadata, PFC.metadata, VTA_SnR.metadata)

write.csv(allcellmetadata, file = "./amb_PT_allcellmetadata.ZW.20240614.csv", sep = ",")

test <- read.csv("./amb_PT_allcellmetadata.ZW.20240614.csv", row.names = 1)

table(test$umap_filter)
table(test$L5r)

test["A01:C6:JM:09",]

table(test$sublib)

table(allcellmetadata$brainregion,allcellmetadata$umap_filter)


HCa.clean <- mcreadRDS(file = "./HCa.clean.rna.rds")
###########################
umap.plot1 <- umap_plots$`06 CTX-CGE GABA`
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, cells_to_highlight)
Idents(VTA_SnR.clean) <- "umap_filter"
Idents(VTA_SnR.clean, cells = select.cells) <- "keep"
VTA_SnR.clean[["umap_filter"]] <- Idents(VTA_SnR.clean)





check <- subset(x = HCa.clean, subset = class_id_label =="03 OB-CR Glut")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")

checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(HCa.clean@meta.data$subclass_id_label)
DimPlot(HCa.clean, group.by = "subclass_id_label", reduction = "L1UMAP")
DimPlot(HCa.clean, group.by = "umap_filter", reduction = "L1UMAP")
DimPlot(HCa.clean, group.by = "class_id_label", reduction = "L1UMAP")



################
umap.plot1 <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, rownames(check@meta.data))






check <- subset(x = HCa.clean, subset = class_id_label =="04 DG-IMN Glut")
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.2, raster = FALSE)
DimPlot(check, group.by = "supertype_id_label", reduction = "L1UMAP", pt.size = 0.2, raster = FALSE)
DimPlot(check, group.by = "umap_filter", reduction = "L1UMAP", pt.size = 0.2, raster = FALSE)

DimPlot(check, group.by = "supertype_id_label", split.by = "umap_filter",reduction = "L1UMAP", pt.size = 0.2, raster = FALSE)




checkplot <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP", pt.size = 0.1)
HoverLocator(plot = checkplot, information = FetchData(check, vars = c("class_id_label", "subclass_id_label","supertype_id_label")))

table(HCa.clean@meta.data$subclass_id_label)
DimPlot(HCa.clean, group.by = "subclass_id_label", reduction = "L1UMAP")
DimPlot(HCa.clean, group.by = "umap_filter", reduction = "L1UMAP")
DimPlot(HCa.clean, group.by = "class_id_label", reduction = "L1UMAP")



################
umap.plot1 <- DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")
select.cells <- CellSelector(plot = umap.plot1)
select.cells <- intersect(select.cells, rownames(check@meta.data))

check
DimPlot(check, group.by = "subclass_id_label", reduction = "L1UMAP")
DimPlot(check, group.by = "supertype_id_label", reduction = "L1UMAP")

###################################
# 20240618 Check no-reference cells from integration

AMY.clean <- mcreadRDS(file = "./AMY.clean.rna.rds")

highlightcells <- read.delim("./un_coembed/AMY.barcode.txt", header = FALSE)
highlightcells <- highlightcells$V1


AMY.clean[["coembed_filter"]] <- "coembed"
Idents(AMY.clean) <- "coembed_filter"
Idents(AMY.clean, cells = highlightcells) <- "uncoembed"
AMY.clean[["coembed_filter"]] <- Idents(AMY.clean)



p1 <- DimPlot(AMY.clean, cells.highlight = highlightcells, reduction = "L1UMAP") + 
  ggtitle("AMY_L1UMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(AMY.clean, cells.highlight = highlightcells, reduction = "umap") + 
  ggtitle("AMY_SeuratUMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- DimPlot(AMY.clean, group.by = "coembed_filter",split.by = "umap_filter", reduction = "L1UMAP", order = "uncoembed", cols = c("coembed" = "lightgrey", "uncoembed" = "red2")) + 
  ggtitle("AMY_L1UMAP_class_id_label")+
  theme(plot.title = element_text(hjust = 0.5))


(p1 | p2) /
  p3



table(AMY.clean@meta.data$umap_filter, AMY.clean@meta.data$coembed_filter)



rm(AMY.clean, highlightcells, p1, p2,p3);gc()


###################################################

CPU.clean <- mcreadRDS(file = "./CPU.clean.rna.rds")
highlightcells <- read.delim("./un_coembed/CPU.barcode.1.txt", header = FALSE)
highlightcells <- highlightcells$V1


CPU.clean[["coembed_filter"]] <- "coembed"
Idents(CPU.clean) <- "coembed_filter"
Idents(CPU.clean, cells = highlightcells) <- "uncoembed"
CPU.clean[["coembed_filter"]] <- Idents(CPU.clean)



p1 <- DimPlot(CPU.clean, cells.highlight = highlightcells, reduction = "L1UMAP") + 
  ggtitle("CPU_L1UMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(CPU.clean, cells.highlight = highlightcells, reduction = "umap") + 
  ggtitle("CPU_SeuratUMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- DimPlot(CPU.clean, group.by = "coembed_filter",split.by = "umap_filter", reduction = "L1UMAP", order = "uncoembed", cols = c("coembed" = "lightgrey", "uncoembed" = "red2")) + 
  ggtitle("CPU_L1UMAP_class_id_label")+
  theme(plot.title = element_text(hjust = 0.5))


(p1 | p2) /
  p3

table(CPU.clean@meta.data$umap_filter, CPU.clean@meta.data$coembed_filter)


highlightcells <- read.delim("./un_coembed/CPU.barcode.2.txt", header = FALSE)
highlightcells <- highlightcells$V1


CPU.clean[["coembed_filter"]] <- "coembed"
Idents(CPU.clean) <- "coembed_filter"
Idents(CPU.clean, cells = highlightcells) <- "uncoembed"
CPU.clean[["coembed_filter"]] <- Idents(CPU.clean)



p1 <- DimPlot(CPU.clean, cells.highlight = highlightcells, reduction = "L1UMAP") + 
  ggtitle("CPU_L1UMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(CPU.clean, cells.highlight = highlightcells, reduction = "umap") + 
  ggtitle("CPU_SeuratUMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- DimPlot(CPU.clean, group.by = "coembed_filter",split.by = "umap_filter", reduction = "L1UMAP", order = "uncoembed", cols = c("coembed" = "lightgrey", "uncoembed" = "red2")) + 
  ggtitle("CPU_L1UMAP_class_id_label")+
  theme(plot.title = element_text(hjust = 0.5))


(p1 | p2) /
  p3

table(CPU.clean@meta.data$umap_filter, CPU.clean@meta.data$coembed_filter)



highlightcells <- read.delim("./un_coembed/CPU.barcode.3.txt", header = FALSE)
highlightcells <- highlightcells$V1


CPU.clean[["coembed_filter"]] <- "coembed"
Idents(CPU.clean) <- "coembed_filter"
Idents(CPU.clean, cells = highlightcells) <- "uncoembed"
CPU.clean[["coembed_filter"]] <- Idents(CPU.clean)



p1 <- DimPlot(CPU.clean, cells.highlight = highlightcells, reduction = "L1UMAP") + 
  ggtitle("CPU_L1UMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(CPU.clean, cells.highlight = highlightcells, reduction = "umap") + 
  ggtitle("CPU_SeuratUMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- DimPlot(CPU.clean, group.by = "coembed_filter",split.by = "umap_filter", reduction = "L1UMAP", order = "uncoembed", cols = c("coembed" = "lightgrey", "uncoembed" = "red2")) + 
  ggtitle("CPU_L1UMAP_class_id_label")+
  theme(plot.title = element_text(hjust = 0.5))


(p1 | p2) /
  p3

table(CPU.clean@meta.data$umap_filter, CPU.clean@meta.data$coembed_filter)

rm(CPU.clean, highlightcells, p1, p2,p3);gc()


############################################


ERC.clean <- mcreadRDS(file = "./ERC.clean.rna.rds")


highlightcells <- read.delim("./un_coembed/ERC.barcode.txt", header = FALSE)
highlightcells <- highlightcells$V1


ERC.clean[["coembed_filter"]] <- "coembed"
Idents(ERC.clean) <- "coembed_filter"
Idents(ERC.clean, cells = highlightcells) <- "uncoembed"
ERC.clean[["coembed_filter"]] <- Idents(ERC.clean)



p1 <- DimPlot(ERC.clean, cells.highlight = highlightcells, reduction = "L1UMAP") + 
  ggtitle("ERC_L1UMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(ERC.clean, cells.highlight = highlightcells, reduction = "umap") + 
  ggtitle("ERC_SeuratUMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- DimPlot(ERC.clean, group.by = "coembed_filter",split.by = "umap_filter", reduction = "L1UMAP", order = "uncoembed", cols = c("coembed" = "lightgrey", "uncoembed" = "red2")) + 
  ggtitle("ERC_L1UMAP_class_id_label")+
  theme(plot.title = element_text(hjust = 0.5))


(p1 | p2) /
  p3



table(ERC.clean@meta.data$umap_filter, ERC.clean@meta.data$coembed_filter)



rm(ERC.clean, highlightcells, p1, p2,p3);gc()



####################################

HCa.clean <- mcreadRDS(file = "./HCa.clean.rna.rds")

highlightcells <- read.delim("./un_coembed/HIP.barcode.1.txt", header = FALSE)
highlightcells <- highlightcells$V1


HCa.clean[["coembed_filter"]] <- "coembed"
Idents(HCa.clean) <- "coembed_filter"
Idents(HCa.clean, cells = highlightcells) <- "uncoembed"
HCa.clean[["coembed_filter"]] <- Idents(HCa.clean)



p1 <- DimPlot(HCa.clean, cells.highlight = highlightcells, reduction = "L1UMAP") + 
  ggtitle("HCa_L1UMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(HCa.clean, cells.highlight = highlightcells, reduction = "umap") + 
  ggtitle("HCa_SeuratUMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- DimPlot(HCa.clean, group.by = "coembed_filter",split.by = "umap_filter", reduction = "L1UMAP", order = "uncoembed", cols = c("coembed" = "lightgrey", "uncoembed" = "red2")) + 
  ggtitle("HCa_L1UMAP_class_id_label")+
  theme(plot.title = element_text(hjust = 0.5))


(p1 | p2) /
  p3



table(HCa.clean@meta.data$umap_filter, HCa.clean@meta.data$coembed_filter)




highlightcells <- read.delim("./un_coembed/HIP.barcode.2.txt", header = FALSE)
highlightcells <- highlightcells$V1


HCa.clean[["coembed_filter"]] <- "coembed"
Idents(HCa.clean) <- "coembed_filter"
Idents(HCa.clean, cells = highlightcells) <- "uncoembed"
HCa.clean[["coembed_filter"]] <- Idents(HCa.clean)



p1 <- DimPlot(HCa.clean, cells.highlight = highlightcells, reduction = "L1UMAP") + 
  ggtitle("HCa_L1UMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(HCa.clean, cells.highlight = highlightcells, reduction = "umap") + 
  ggtitle("HCa_SeuratUMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- DimPlot(HCa.clean, group.by = "coembed_filter",split.by = "umap_filter", reduction = "L1UMAP", order = "uncoembed", cols = c("coembed" = "lightgrey", "uncoembed" = "red2")) + 
  ggtitle("HCa_L1UMAP_class_id_label")+
  theme(plot.title = element_text(hjust = 0.5))


(p1 | p2) /
  p3



table(HCa.clean@meta.data$umap_filter, HCa.clean@meta.data$coembed_filter)





highlightcells <- read.delim("./un_coembed/HIP.barcode.3.txt", header = FALSE)
highlightcells <- highlightcells$V1


HCa.clean[["coembed_filter"]] <- "coembed"
Idents(HCa.clean) <- "coembed_filter"
Idents(HCa.clean, cells = highlightcells) <- "uncoembed"
HCa.clean[["coembed_filter"]] <- Idents(HCa.clean)



p1 <- DimPlot(HCa.clean, cells.highlight = highlightcells, reduction = "L1UMAP") + 
  ggtitle("HCa_L1UMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(HCa.clean, cells.highlight = highlightcells, reduction = "umap") + 
  ggtitle("HCa_SeuratUMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- DimPlot(HCa.clean, group.by = "coembed_filter",split.by = "umap_filter", reduction = "L1UMAP", order = "uncoembed", cols = c("coembed" = "lightgrey", "uncoembed" = "red2")) + 
  ggtitle("HCa_L1UMAP_class_id_label")+
  theme(plot.title = element_text(hjust = 0.5))


(p1 | p2) /
  p3



table(HCa.clean@meta.data$umap_filter, HCa.clean@meta.data$coembed_filter)







rm(HCa.clean, highlightcells, p1, p2,p3);gc()
###################################


HCp.clean <- mcreadRDS(file = "./HCp.clean.rna.rds")


highlightcells <- read.delim("./un_coembed/HIP.barcode.1.txt", header = FALSE)
highlightcells <- highlightcells$V1


HCp.clean[["coembed_filter"]] <- "coembed"
Idents(HCp.clean) <- "coembed_filter"
Idents(HCp.clean, cells = highlightcells) <- "uncoembed"
HCp.clean[["coembed_filter"]] <- Idents(HCp.clean)



p1 <- DimPlot(HCp.clean, cells.highlight = highlightcells, reduction = "L1UMAP") + 
  ggtitle("HCp_L1UMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(HCp.clean, cells.highlight = highlightcells, reduction = "umap") + 
  ggtitle("HCp_SeuratUMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- DimPlot(HCp.clean, group.by = "coembed_filter",split.by = "umap_filter", reduction = "L1UMAP", order = "uncoembed", cols = c("coembed" = "lightgrey", "uncoembed" = "red2")) + 
  ggtitle("HCp_L1UMAP_class_id_label")+
  theme(plot.title = element_text(hjust = 0.5))


(p1 | p2) /
  p3



table(HCp.clean@meta.data$umap_filter, HCp.clean@meta.data$coembed_filter)




highlightcells <- read.delim("./un_coembed/HIP.barcode.2.txt", header = FALSE)
highlightcells <- highlightcells$V1


HCp.clean[["coembed_filter"]] <- "coembed"
Idents(HCp.clean) <- "coembed_filter"
Idents(HCp.clean, cells = highlightcells) <- "uncoembed"
HCp.clean[["coembed_filter"]] <- Idents(HCp.clean)



p1 <- DimPlot(HCp.clean, cells.highlight = highlightcells, reduction = "L1UMAP") + 
  ggtitle("HCp_L1UMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(HCp.clean, cells.highlight = highlightcells, reduction = "umap") + 
  ggtitle("HCp_SeuratUMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- DimPlot(HCp.clean, group.by = "coembed_filter",split.by = "umap_filter", reduction = "L1UMAP", order = "uncoembed", cols = c("coembed" = "lightgrey", "uncoembed" = "red2")) + 
  ggtitle("HCp_L1UMAP_class_id_label")+
  theme(plot.title = element_text(hjust = 0.5))


(p1 | p2) /
  p3



table(HCp.clean@meta.data$umap_filter, HCp.clean@meta.data$coembed_filter)


highlightcells <- read.delim("./un_coembed/HIP.barcode.3.txt", header = FALSE)
highlightcells <- highlightcells$V1


HCp.clean[["coembed_filter"]] <- "coembed"
Idents(HCp.clean) <- "coembed_filter"
Idents(HCp.clean, cells = highlightcells) <- "uncoembed"
HCp.clean[["coembed_filter"]] <- Idents(HCp.clean)



p1 <- DimPlot(HCp.clean, cells.highlight = highlightcells, reduction = "L1UMAP") + 
  ggtitle("HCp_L1UMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(HCp.clean, cells.highlight = highlightcells, reduction = "umap") + 
  ggtitle("HCp_SeuratUMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- DimPlot(HCp.clean, group.by = "coembed_filter",split.by = "umap_filter", reduction = "L1UMAP", order = "uncoembed", cols = c("coembed" = "lightgrey", "uncoembed" = "red2")) + 
  ggtitle("HCp_L1UMAP_class_id_label")+
  theme(plot.title = element_text(hjust = 0.5))


(p1 | p2) /
  p3



table(HCp.clean@meta.data$umap_filter, HCp.clean@meta.data$coembed_filter)




rm(HCp.clean, highlightcells, p1, p2,p3);gc()






###################################
HYP.clean <- mcreadRDS(file = "./HYP.clean.rna.rds")

highlightcells <- read.delim("./un_coembed/HYP.barcode.txt", header = FALSE)
highlightcells <- highlightcells$V1


HYP.clean[["coembed_filter"]] <- "coembed"
Idents(HYP.clean) <- "coembed_filter"
Idents(HYP.clean, cells = highlightcells) <- "uncoembed"
HYP.clean[["coembed_filter"]] <- Idents(HYP.clean)



p1 <- DimPlot(HYP.clean, cells.highlight = highlightcells, reduction = "L1UMAP") + 
  ggtitle("HYP_L1UMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(HYP.clean, cells.highlight = highlightcells, reduction = "umap") + 
  ggtitle("HYP_SeuratUMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- DimPlot(HYP.clean, group.by = "coembed_filter",split.by = "umap_filter", reduction = "L1UMAP", order = "uncoembed", cols = c("coembed" = "lightgrey", "uncoembed" = "red2")) + 
  ggtitle("HYP_L1UMAP_class_id_label")+
  theme(plot.title = element_text(hjust = 0.5))


(p1 | p2) /
  p3



table(HYP.clean@meta.data$umap_filter, HYP.clean@meta.data$coembed_filter)



rm(HYP.clean, highlightcells, p1, p2,p3);gc()







################################
NAC.clean <- mcreadRDS(file = "./NAC.clean.rna.rds")

highlightcells <- read.delim("./un_coembed/NAC.barcode.txt", header = FALSE)
highlightcells <- highlightcells$V1


NAC.clean[["coembed_filter"]] <- "coembed"
Idents(NAC.clean) <- "coembed_filter"
Idents(NAC.clean, cells = highlightcells) <- "uncoembed"
NAC.clean[["coembed_filter"]] <- Idents(NAC.clean)



p1 <- DimPlot(NAC.clean, cells.highlight = highlightcells, reduction = "L1UMAP") + 
  ggtitle("NAC_L1UMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(NAC.clean, cells.highlight = highlightcells, reduction = "umap") + 
  ggtitle("NAC_SeuratUMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- DimPlot(NAC.clean, group.by = "coembed_filter",split.by = "umap_filter", reduction = "L1UMAP", order = "uncoembed", cols = c("coembed" = "lightgrey", "uncoembed" = "red2")) + 
  ggtitle("NAC_L1UMAP_class_id_label")+
  theme(plot.title = element_text(hjust = 0.5))


(p1 | p2) /
  p3



table(NAC.clean@meta.data$umap_filter, NAC.clean@meta.data$coembed_filter)



rm(NAC.clean, highlightcells, p1, p2,p3);gc()



#################################
PFC.clean <- mcreadRDS(file = "./PFC.clean.rna.rds")

highlightcells <- read.delim("./un_coembed/PFC.barcode.1.txt", header = FALSE)
highlightcells <- highlightcells$V1


PFC.clean[["coembed_filter"]] <- "coembed"
Idents(PFC.clean) <- "coembed_filter"
Idents(PFC.clean, cells = highlightcells) <- "uncoembed"
PFC.clean[["coembed_filter"]] <- Idents(PFC.clean)



p1 <- DimPlot(PFC.clean, cells.highlight = highlightcells, reduction = "L1UMAP") + 
  ggtitle("PFC_L1UMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(PFC.clean, cells.highlight = highlightcells, reduction = "umap") + 
  ggtitle("PFC_SeuratUMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- DimPlot(PFC.clean, group.by = "coembed_filter",split.by = "umap_filter", reduction = "L1UMAP", order = "uncoembed", cols = c("coembed" = "lightgrey", "uncoembed" = "red2")) + 
  ggtitle("PFC_L1UMAP_class_id_label")+
  theme(plot.title = element_text(hjust = 0.5))


(p1 | p2) /
  p3



table(PFC.clean@meta.data$umap_filter, PFC.clean@meta.data$coembed_filter)


highlightcells <- read.delim("./un_coembed/PFC.barcode.2.txt", header = FALSE)
highlightcells <- highlightcells$V1


PFC.clean[["coembed_filter"]] <- "coembed"
Idents(PFC.clean) <- "coembed_filter"
Idents(PFC.clean, cells = highlightcells) <- "uncoembed"
PFC.clean[["coembed_filter"]] <- Idents(PFC.clean)



p1 <- DimPlot(PFC.clean, cells.highlight = highlightcells, reduction = "L1UMAP") + 
  ggtitle("PFC_L1UMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(PFC.clean, cells.highlight = highlightcells, reduction = "umap") + 
  ggtitle("PFC_SeuratUMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- DimPlot(PFC.clean, group.by = "coembed_filter",split.by = "umap_filter", reduction = "L1UMAP", order = "uncoembed", cols = c("coembed" = "lightgrey", "uncoembed" = "red2")) + 
  ggtitle("PFC_L1UMAP_class_id_label")+
  theme(plot.title = element_text(hjust = 0.5))


(p1 | p2) /
  p3



table(PFC.clean@meta.data$umap_filter, PFC.clean@meta.data$coembed_filter)






highlightcells <- read.delim("./un_coembed/PFC.barcode.3.txt", header = FALSE)
highlightcells <- highlightcells$V1


PFC.clean[["coembed_filter"]] <- "coembed"
Idents(PFC.clean) <- "coembed_filter"
Idents(PFC.clean, cells = highlightcells) <- "uncoembed"
PFC.clean[["coembed_filter"]] <- Idents(PFC.clean)



p1 <- DimPlot(PFC.clean, cells.highlight = highlightcells, reduction = "L1UMAP") + 
  ggtitle("PFC_L1UMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(PFC.clean, cells.highlight = highlightcells, reduction = "umap") + 
  ggtitle("PFC_SeuratUMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- DimPlot(PFC.clean, group.by = "coembed_filter",split.by = "umap_filter", reduction = "L1UMAP", order = "uncoembed", cols = c("coembed" = "lightgrey", "uncoembed" = "red2")) + 
  ggtitle("PFC_L1UMAP_class_id_label")+
  theme(plot.title = element_text(hjust = 0.5))


(p1 | p2) /
  p3



table(PFC.clean@meta.data$umap_filter, PFC.clean@meta.data$coembed_filter)


rm(PFC.clean, highlightcells, p1, p2,p3);gc()

############################
VTA_SnR.clean <- mcreadRDS(file = "./VTA_SnR.clean.rna.rds")
highlightcells <- read.delim("./un_coembed/VTA.barcode.1.txt", header = FALSE)
highlightcells <- highlightcells$V1


VTA_SnR.clean[["coembed_filter"]] <- "coembed"
Idents(VTA_SnR.clean) <- "coembed_filter"
Idents(VTA_SnR.clean, cells = highlightcells) <- "uncoembed"
VTA_SnR.clean[["coembed_filter"]] <- Idents(VTA_SnR.clean)



p1 <- DimPlot(VTA_SnR.clean, cells.highlight = highlightcells, reduction = "L1UMAP") + 
  ggtitle("VTA_SnR_L1UMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(VTA_SnR.clean, cells.highlight = highlightcells, reduction = "umap") + 
  ggtitle("VTA_SnR_SeuratUMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- DimPlot(VTA_SnR.clean, group.by = "coembed_filter",split.by = "umap_filter", reduction = "L1UMAP", order = "uncoembed", cols = c("coembed" = "lightgrey", "uncoembed" = "red2")) + 
  ggtitle("VTA_SnR_L1UMAP_class_id_label")+
  theme(plot.title = element_text(hjust = 0.5))


(p1 | p2) /
  p3



table(VTA_SnR.clean@meta.data$umap_filter, VTA_SnR.clean@meta.data$coembed_filter)




highlightcells <- read.delim("./un_coembed/VTA.barcode.2.txt", header = FALSE)
highlightcells <- highlightcells$V1


VTA_SnR.clean[["coembed_filter"]] <- "coembed"
Idents(VTA_SnR.clean) <- "coembed_filter"
Idents(VTA_SnR.clean, cells = highlightcells) <- "uncoembed"
VTA_SnR.clean[["coembed_filter"]] <- Idents(VTA_SnR.clean)



p1 <- DimPlot(VTA_SnR.clean, cells.highlight = highlightcells, reduction = "L1UMAP") + 
  ggtitle("VTA_SnR_L1UMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(VTA_SnR.clean, cells.highlight = highlightcells, reduction = "umap") + 
  ggtitle("VTA_SnR_SeuratUMAP_uncoembedded")+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- DimPlot(VTA_SnR.clean, group.by = "coembed_filter",split.by = "umap_filter", reduction = "L1UMAP", order = "uncoembed", cols = c("coembed" = "lightgrey", "uncoembed" = "red2")) + 
  ggtitle("VTA_SnR_L1UMAP_class_id_label")+
  theme(plot.title = element_text(hjust = 0.5))


(p1 | p2) /
  p3



table(VTA_SnR.clean@meta.data$umap_filter, VTA_SnR.clean@meta.data$coembed_filter)



rm(VTA_SnR.clean, highlightcells, p1, p2,p3);gc()
