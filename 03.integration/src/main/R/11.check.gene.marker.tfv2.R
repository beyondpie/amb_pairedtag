library(tidyverse)
library(Seurat)
Sys.setenv("_R_USE_PIPEBIND_" = "true")

# * meta
projd <- "/mnt/tscc2/szu/projects/amb_pairedtag"
setwd(projd)
workd <- file.path(projd, "03.integration")
ptMetaCheckv2 <- file.path(workd, "out", "tfneu_vf_region_cca_k5", 
                           "amb_PT_allcellmetadata.ZW.20240614.csv") |>
  data.table::fread(file = _, sep = ",", 
                    data.table = FALSE, header = TRUE)
final_barcodes <- with(ptMetaCheckv2, V1[umap_filter == "keep"])
rownames(ptMetaCheckv2) <- ptMetaCheckv2$V1

outfigd <- file.path(workd, "out", "tfneu_vf_region_cca_k5")

# table(ptMetaCheckv2$umap_filter)
# table(with(ptMetaCheckv2, umap_filter[brainregion == "AMY"]))

# * get gene expression for paired-tag data on supertype level
#   and subclass level
ptseud <- file.path(projd, "data", "pairedtag_seurat")
nn_ptseu <- readRDS(file.path(ptseud, "ptRNA.nn.k8.L5ds30.rds"))
neu_ptseu <- readRDS(file.path(ptseud, "ptRNA.neu.k8.L5ds50.rds"))
# merge two seurats
ptseu <- merge(nn_ptseu, neu_ptseu,
               add.cell.ids = NULL)
barcodes_ds <- colnames(ptseu) |>
  x => x[x %in% final_barcodes]
ptseu <- ptseu[ , colnames(ptseu) %in% final_barcodes]
ptseu$supertype_id_label <- ptMetaCheckv2[colnames(ptseu), "supertype_id_label"]
ptseu$subclass_id_label <- ptMetaCheckv2[colnames(ptseu), "subclass_id_label"]

pseudo_sp <- AggregateExpression(ptseu, group.by = "supertype_id_label", 
                                 return.seurat = TRUE)
pseudo_sc <- AggregateExpression(ptseu, group.by = "subclass_id_label",
                                 return.seurat = TRUE)

get_id <- function(renamed_strs) {
  vapply(renamed_strs, \(s) {
    str_split_1(s, " ")[1] |>
      gsub("g", "", x = _) |>
      as.integer()
  }, 1L)
}
spid_pseudo <- get_id(pseudo_sp$supertype_id_label)
scid_pseudo <- get_id(pseudo_sc$subclass_id_label)

spid_allen <- get_id(allenclMeta$supertype_id_label)
scid_allen <- get_id(allenclMeta$subclass_id_label)
spid2sp_allen <- data.frame(
  id = paste0("id-", spid_allen),
  sp = allenclMeta$supertype_id_label
) |>
  unique() |>
  x => `rownames<-`(x, x$id)
scid2sc_allen <- data.frame(
  id = paste0("id-", scid_allen),
  sc = allenclMeta$subclass_id_label
) |>
  unique() |>
  x => `rownames<-`(x, x$id)

pseudo_sp$sp <- spid2sp_allen[
  paste0("id-", spid_pseudo),
  "sp"
]
pseudo_sc$sc <- scid2sc_allen[
  paste0("id-", scid_pseudo),
  "sc"
]
# order based on name id
pseudo_sp <- pseudo_sp[ , order(pseudo_sp$sp)]
pseudo_sc <- pseudo_sc[ , order(pseudo_sc$sc)]

# do normalize and scale data on ptseu
ptseu <- NormalizeData(ptseu)
ptseu <- ScaleData(ptseu)

a <- ptseu[, c(2,1)]

# * select marker genes
allenclMeta <- file.path(projd, "meta",
                         "AIT21_annotation_freeze_081523.tsv") |>
  data.table::fread(file = _, sep = "\t", header = TRUE, 
                    data.table = FALSE)

get_allen_top_marker_by_group <- function(onAllenLevel, top = 3) {
  rank_markers <- function(mgs) {
    g <- str_split_1(mgs, ",")
    table(g) |>
      as.data.frame(x = _, stringsAsFactors = FALSE) |>
      x => `rownames<-`(x, x$genes) |>
      x => x[order(x$Freq, decreasing = TRUE), ]
  }
  
  get_topk_markers <- function(mgs, top = 3) {
    rank_markers(mgs) |>
      head(x = _, n = top)
  }
  
  sc2mgs <- allenclMeta |>
    group_by(.data[[onAllenLevel]]) |>
    summarise(marker = paste(cluster.markers, collapse = ",")) |>
    as.data.frame()
  
  r <- lapply(seq_len(nrow(sc2mgs)), function(i) {
    get_topk_markers(sc2mgs[i, "marker"], top = top)
  }) |>
    setNames(object = _, nm = sc2mgs[[onAllenLevel]])
  return(r)
}

sp_mg <- get_allen_top_marker_by_group("supertype_id_label", top = 3)
sc_mg <- get_allen_top_marker_by_group("subclass_id_label", top = 3)

# common markers
commonGenes <- data.table::fread(
  file = file.path(
    projd, "meta",
    "common_gene_markers.csv"
  ),
  sep = ",", header = TRUE, data.table = FALSE
)

slt_cmn_gene <- commonGenes[1:13, "gene"] |>
  x => x[x %in% rownames(pseudo_sc)]

slt_sp_gene <- lapply(pseudo_sp$sp, \(sp) {
  sp_mg[[sp]][, "g"]
}) |> do.call("c", args = _) |>
  x => x[x %in% rownames(pseudo_sp)]

slt_sc_gene <- lapply(pseudo_sc$sc, \(sc) {
  sc_mg[[sc]][ , "g"]
}) |> do.call("c", args = _) |>
  x => x[x %in% rownames(pseudo_sc)] |>
  unique()

# * plot
mytheme <- theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  plot.title = element_text(colour = "black", hjust = 0.5,
                            size = 15),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_text(colour = "black", size = 12),
  axis.line = element_line(colour = "black"),
  legend.position = "right",
  legend.title = element_text(colour = "black", size = 13),
  legend.text = element_text(colour = "black", size = 12),
)

# * plot gene expression on paired-tag pseudo-bulk level
p_pseudo_sc <- DotPlot(pseudo_sc, features = unique(c(slt_cmn_gene, slt_sc_gene)),
                cols = c("white", "blue"),
                dot.scale = 0.5, col.min = 0)
p_pseudo_sc + mytheme

p_pseudo_sp <- DotPlot(pseudo_sp, features = unique(c(slt_cmn_gene, slt_sp_gene)),
                       cols = c("white", "blue"),
                       dot.scale = 0.5, col.min = 0)
p_pseudo_sp + mytheme

# * plot gene expression on paired-tag single-cell level
p_singlecell_sp <- DotPlot(ptseu, 
                           group.by = "supertype_id_label",
                           features = unique(c(slt_cmn_gene, slt_sp_gene)),
                           cols = c("white", "blue"),
                           dot.scale = 2, col.min = 0) +
  ggtitle("PairedTag Single-cell RNA-seq Gene Expression Grouped by SuperType") + 
  mytheme

ggsave(plot = p_singlecell_sp,
       filename = file.path(outfigd, "pt.singlecell.supertype.marker.pdf"),
       width = 25, height = 10)

p_singlecell_sc <- DotPlot(ptseu, 
                           group.by = "subclass_id_label",
                           features = unique(c(slt_cmn_gene, slt_sc_gene)),
                           cols = c("white", "blue"),
                           dot.scale = 2, col.min = 0) + 
  ggtitle("PairedTag Single-cell RNA-seq Gene Expression Grouped by Subclass") + 
  mytheme
ggsave(plot = p_singlecell_sc, 
       filename = file.path(outfigd, "pt.singlecell.subclass.marker.pdf"),
       width = 15, height = 10)

# * plot gene expression on Allen's single-cell level
allen_nn_seu <- file.path(projd, "data", "allen_seurat",
                          "allen.10xv3.pt.regions.nn.imn.k8.cl.ds1000.rds") |>
  readRDS(file = _)
allen_neu_seu <- file.path(projd, "data", "allen_seurat",
                           "allen.10xv3.pt.regions.neu.imn.k8.cl.ds100.rds") |>
  readRDS(file = _)

allen_neu_seu <- allen_neu_seu[, 
                               !(colnames(allen_neu_seu) %in% colnames(allen_nn_seu))]

allen_seu <- merge(allen_nn_seu, allen_neu_seu,  add.cell.ids = NULL)
allen_seu <- NormalizeData(allen_seu) |>
  ScaleData(object = _ )
allen_seu <- subset(allen_seu, supertype_id_label %in% ptseu$supertype_id_label)
p_allen_sp <- DotPlot(allen_seu, 
                           group.by = "supertype_id_label",
                           features = unique(c(slt_cmn_gene, slt_sp_gene)),
                           cols = c("white", "blue"),
                           dot.scale = 2, col.min = 0) +
  ggtitle("Allen Single-cell RNA-seq Gene Expression Grouped by SuperType") + 
  mytheme

ggsave(plot = p_allen_sp,
       filename = file.path(outfigd, "allen.singlecell.supertype.marker.pdf"),
       width = 25, height = 10)

p_allen_sc <- DotPlot(allen_seu, 
                           group.by = "subclass_id_label",
                           features = unique(c(slt_cmn_gene, slt_sc_gene)),
                           cols = c("white", "blue"),
                           dot.scale = 2, col.min = 0) + 
  ggtitle("Allen Single-cell RNA-seq Gene Expression Grouped by Subclass") + 
  mytheme
ggsave(plot = p_allen_sc, 
       filename = file.path(outfigd, "allen.singlecell.subclass.marker.pdf"),
       width = 15, height = 10)

