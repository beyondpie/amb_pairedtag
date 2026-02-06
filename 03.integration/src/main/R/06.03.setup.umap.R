library(tidyverse)
library(Seurat)
library(ShinyCell)
library(rsconnect)


# shared meta
projd <- here::here()
clprefix <- "cl-"
onAllenLevel <- "subclass_id_label"

# neuron seurat
tfSumdir <- file.path(projd, "03.integration", "out",
  "tfneu_k8_cca_k5")
tfsumfnm <- file.path(tfSumdir, "sum_tfneu_k8_cca_k5.seurat.rds")
cnssfnm <- file.path(tfSumdir, "sum_tfneu.consensus.mat.rds")
ptseu <- readRDS(file.path(projd, "data/pairedtag_seurat",
  "ptRNA.neu.k8.L5ds50.rds"))

# 0. load tf result
tfsum <- readRDS(tfsumfnm)
tfgenes <- rownames(tfsum)

cnss <- readRDS(cnssfnm)
tfL5s <- colnames(cnss$sc)
tfscs <- if (grepl("supertype", onAllenLevel)) {
  rownames(cnss$sp)
} else {
  rownames(cnss$sc)
}

# * allen cluster meta
allenclMetafnm <- file.path(projd, "meta",
  "AIT21_annotation_freeze_081523.tsv")
allenclMeta <- data.table::fread(
  file = allenclMetafnm, sep = "\t", header = TRUE,
  data.table = FALSE
)
rownames(allenclMeta) <- paste0(clprefix, allenclMeta$cl)
# filter by subclass we focus on
allenclMeta <- allenclMeta[allenclMeta[[onAllenLevel]] %in% tfscs, ]

# paired tag meta data
ptcellMeta <- file.path(
  projd, "meta",
  "pairedtag.cell.meta.all.csv"
) |>
  data.table::fread(
    file = _, sep = ",", header = TRUE,
    data.table = FALSE
  ) |>
  x => `rownames<-`(x, x$barcode) |>
  x => x[paste0("L5-", x$L1_2_3_4_5) %in% tfL5s, ]

# merge subtype, subclass info
# add UMAP
ptseu <- NormalizeData(ptseu)
ptseu <- FindVariableFeatures(ptseu, nfeatures = 3000)
ptseu <- ScaleData(ptseu, features = rownames(ptseu))
ptseu <- RunPCA(ptseu, features = VariableFeatures(ptseu),
  npcs = 50)
ptseu <- RunUMAP(ptseu, dims = 1:50)
# add supertype info
L5_sp <- tfsum@meta.data[!tfsum$isRef, ] |>
  group_by(L5) |>
  summarise(sp = unique(supertype_id_label)) |>
  as.data.frame() |>
  x => `rownames<-`(x, paste0("L5-", x$L5))
ptcellMeta$supertype_id_label <- L5_sp[
  paste0("L5-", ptcellMeta$L1_2_3_4_5), "sp"]
ptcellMeta$cl <- L5_cl[
  paste0("L5-", ptcellMeta$L1_2_3_4_5), "cl"]

a <- factor(ptcellMeta[colnames(ptseu), "supertype_id_label"])
b <- factor(ptcellMeta[colnames(ptseu), "modality"])
ptseu <- AddMetaData(object = ptseu,
  a, col.name = "cluster")
ptseu <- AddMetaData(object = ptseu,
  b, col.name = "DNAmodality")



# * ShinyCell
scConf <- createConfig(ptseu,
  meta.to.include = c("DNAmodality", "cluster")
)
makeShinyApp(
  ptseu, scConf, gene.mapping = TRUE,
  shiny.title = "Pairted-Tag Neuronal Cells"
)

rsconnect::setAccountInfo(name='beyondpie',
			  token='35EC037F341E503ABC1CAD7B2F4E54F9',
			  secret='JuQTMIGxI3E3UdqQNFz4jTtI4bFkzdlwmvRlF79J')
rsconnect::deployApp(appDir = 'shinyApp',
  appName = "Paired-Tag_Neuron_Seurat",
  launch.browser = FALSE)
