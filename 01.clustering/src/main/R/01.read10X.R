library(tidyverse)
library(Seurat)
library(Matrix)

projdir <- here::here()
metadir <- file.path(projdir, "meta")
workdir <- file.path(projdir, "01.clustering")
rscdir <- file.path(workdir, "src/main/resource")
outdir <- file.path(workdir, "out")

# * main
# load seurat files
scRNAFileRecord <- file.path(rscdir, "scRNAseq_samples.passQC.231120.csv")

scRNAFiles <- data.table::fread(
  file = scRNAFileRecord, header = TRUE, sep = ",", data.table = FALSE)
seuList <- lapply(seq_len(nrow(scRNAFiles)), \(i) {
  message("Current exp: ", scRNAFiles$exp[i])
  f <- scRNAFiles$seurat_new[i]
  message("Load Seurat object: ", f)
  readRDS(f)
})
names(seuList) <- scRNAFiles$exp

# merge meta data
metaList <- lapply(seq_along(seuList), \(i) {
  ## expmnt <- names(seuList)[i]
  m <- seuList[[i]]@meta.data
  m$barcode <- rownames(m)
  ## rownames(m) <- paste(expmnt, rownames(m), sep = "_")
  return(m)
})
names(metaList) <- names(seuList)

# * fix fields in Exp13 and Exp14
# ** Exp13
# "replicate" should be the metadata slot for "rep"
metaList$Exp13$rep <- metaList$Exp13$replicate
metaList$Exp13$replicate <- NULL

# ** Exp14
# 1. correction: (oriBarcode -> brainregion): 23-> ERC; 24->AMY
tm <- metaList$Exp14
metaList$Exp14$brainregion[tm$oriBarcode == "23"] <- "ERC"
metaList$Exp14$brainregion[tm$oriBarcode == "24"] <- "AMY"
# with(metaList$Exp14, table(brainregion, oriBarcode))

# 2. rep info below, and replicate now is ok
## replicate <- rep("MaleA",dim(bc)[1])
## replicate[bc[,4]=="02"|bc[,4]=="14"]<-"MaleB"
## replicate[bc[,4]=="11"|bc[,4]=="21"]<-"FemaleA"
## replicate[bc[,4]=="12"|bc[,4]=="22"]<-"FemaleB"
metaList$Exp14$rep <- metaList$Exp14$replicate
metaList$Exp14$replicate <- NULL
# 3. add sex info
sex <- metaList$Exp14$rep |> gsub("[A|B]$", "", x = _)
metaList$Exp14$sex <- sex

saveRDS(metaList, file.path(rscdir, "scRNA.metaList.v1.rds"))

# * merge meta data
metaAll <- do.call(rbind, metaList)
metav1 <- metaAll |>
  dplyr::select(.data = _, -orig.ident)

data.table::fwrite(
  x = metav1,
  file = file.path(projdir, "meta", "pairedtag.meta.v1.csv"),
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


# * merge seurat objects
metaList2 <- lapply(seq_along(metaList), \(i) {
  m <- metaList[[i]]
  m$barcode <- NULL
  return(m)
})
names(metaList2) <- names(metaList)
seuList$Exp13@meta.data <- metaList2$Exp13
seuList$Exp14@meta.data <- metaList2$Exp14

# * merge seurat files
## take too much time and failed
## vector::reserve
## seuAll <- scCustomize::Merge_Seurat_List(
##   list_seurat = seuList,
##   add.cell.ids = names(seuList),
##   merge.data = TRUE,
##   project = "amb_pairedtag_RNA"
## )

countList <- lapply(seq_along(seuList), \(i) {
  seu <- seuList[[i]]
  seu@assays$RNA@counts
})
names(countList) <- names(seuList)
# do.call on this will run forever
# manually merge, fast but have the same error at last Error: vector::reserve
# different experiments have different number of genes
## countAll <- countList[[1]]
## for (i in 2:length(countList)) {
##   message("Merge count from: ", names(countList)[i])
##   countAll <- SeuratObject::RowMergeSparseMatrices(countAll, countList[[i]])
## }

# * check number of genes per seurat
ngene <- vapply(seuList, nrow, 1)

# * select the genes across all the experiment.
genes <- rownames(seuList[[1]])
for(i in 2:length(seuList)) {
  genes <- intersect(genes, rownames(seuList[[i]]))
}

cntv2List <- lapply(seq_along(seuList), \(i) {
  seuList[[i]]@assays$RNA@counts[genes, ]
})
names(cntv2List) <- names(seuList)

# Error in cbind.Matrix(x, y, deparse.level = 0L) : 
#  p[length(p)] cannot exceed 2^31-1
## cntv2All <- cntv2List[[1]]
## for (i in 2:length(cntv2List)) {
##   message("Merge count from: ", names(cntv2List)[i])
#   cntv2All <- SeuratObject::RowMergeSparseMatrices(
##     cntv2All, cntv2List[[i]])
## }

# * merge seurat files v2
# now use 8k genes from Allen's data as variable features.
gene_k8 <- data.table::fread(
  file = file.path(metadir, "AIT21_k8_markers.txt"), header = FALSE,
  data.table = FALSE)$V1
# 7944
var_genes <- intersect(gene_k8, genes)
cntv3List <- lapply(seq_along(seuList), \(i) {
  seuList[[i]]@assays$RNA@counts[var_genes, ]
})
names(cntv3List) <- names(seuList)
cntv3All <- cntv3List[[1]]
for (i in 2:length(cntv3List)) {
  message("Merge count from: ", names(cntv3List)[i])
  cntv3All <- SeuratObject::RowMergeSparseMatrices(
    cntv3All, cntv3List[[i]])
}

s5 <- Seurat::CreateSeuratObject(counts = cntv3All, assay = "RNA")

metaAll2 <- do.call(rbind, metaList2)
rownm_meta <- gsub("^Exp\\d+.", "", rownames(metaAll2))
rownames(metaAll2) <- rownm_meta
metaAll2 <- metaAll2[colnames(s5), ]
s5 <- Seurat::AddMetaData(s5, metadata = metaAll2)
saveRDS(s5,
  file = file.path(outdir, "scRNAseq_k8_all.seurat.rds"))

# * check NA when for the count data
any(is.na(cntv3All))
any(is.infinite(cntv3All))
any(is.nan(cntv3All))
