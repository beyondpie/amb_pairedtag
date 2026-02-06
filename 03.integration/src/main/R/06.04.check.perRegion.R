library(tidyverse)
library(Seurat)
library(ComplexHeatmap)

# shared meta
projd <- here::here()
clprefix <- "cl-"

# * allen region to ours
ptregion2AllenData <- data.table::fread(
  file = file.path(projd, "meta",
    "Allen_reference_data_locator_ZW_20240131.csv"),
  header = TRUE, sep = ",", data.table = FALSE
)
ptregion <- ptregion2AllenData[["Paired-Tag Brain Regions"]]
allen10xv3fnms <- ptregion2AllenData[, 5]
pt2aRegion <- lapply(seq_len(nrow(ptregion2AllenData)), \(i) {
  ptr <- ptregion2AllenData[i, 1]
  tmp <- ptregion2AllenData[i, 5]
  tmp2 <- lapply(str_split_1(tmp, ";"), \(f) {
    gsub("allen_10xv3_|.h5ad| ", "", f) |>
      gsub("__|-", "_", x = _)
  }) |> unlist()
  data.frame(pt = ptr, allen = tmp2)
}) |> do.call(rbind, args = _) |>
  x => unique(x)
ptregions <- unique(pt2aRegion$pt)
aregions <- unique(pt2aRegion$allen)
data.table::fwrite(x = pt2aRegion,
  file = file.path(projd, "meta", "pairedtag2allen.region.csv"),
  sep = ",", col.names = TRUE, row.names = FALSE)

a2pt <- lapply(aregions, \(a) {
  data.frame(
    allen = a,
    pt = with(pt2aRegion, paste(pt[allen %in% a], collapse = ";"))
  )
}) |> do.call(what = "rbind", args = _) |>
  x => `rownames<-`(x, x$allen)


# * allen cellMeta enhancement
# add paired-tag region info
allencellMeta <- data.table::fread(
  file = file.path(
    projd, "data", "allen",
    "allen.10xv3.cell.meta.csv"
  ), sep = ",", header = TRUE,
  data.table = FALSE) |>
  x => `rownames<-`(x, x$barcode)

aptr <- vapply(allencellMeta$roi, \(a) {
  aa <- gsub(" |  ", "", x = a) |>
    gsub("-", "_", x = _)
  if (aa %in% rownames(a2pt)) {
    a2pt[aa, "pt"]
  } else {
    "none"
  }
}, "CPU")
allencellMeta$ptr <- aptr
data.table::fwrite(
  x = allencellMeta,
  file = file.path(projd, "meta", "allen10xv3.meta.with.ptregion.csv"),
  sep = ",",
  col.names = TRUE, row.names = FALSE
)

# * allen supertype and subclass top gene markers
allenclMeta <- data.table::fread(
  file = file.path(projd, "meta", "AIT21_annotation_freeze_081523.tsv"),
  sep = "\t", header = TRUE, data.table = FALSE
) |>
  x => `rownames<-`(x, paste0(clprefix, x$cl))
get_rank_markers <- function(onAllenLevel = "supertype_id_label") {
  rank_markers <- function(mgs) {
    g <- str_split_1(mgs, ",")
    table(g) |>
      as.data.frame(x = _, stringsAsFactors = FALSE) |>
      x => `rownames<-`(x, x$genes) |>
      x => x[order(x$Freq, decreasing = TRUE), ]
  }
  sc2mgs <- allenclMeta |>
    group_by(.data[[onAllenLevel]]) |>
    summarise(marker = paste(cluster.markers, collapse = ",")) |>
    as.data.frame()

  lapply(seq_len(nrow(sc2mgs)), function(i) {
    rank_markers(sc2mgs[i, "marker"])
  }) |>
    setNames(object = _, nm = sc2mgs[[onAllenLevel]])
}
sp2rankmarkers <- get_rank_markers(
  onAllenLevel = "supertype_id_label")
sc2rankmarkers <- get_rank_markers(
  onAllenLevel = "subclass_id_label"
)
saveRDS(
  object = list(sp = sp2rankmarkers, sc = sc2rankmarkers),
  file = file.path(projd, "meta", "allen_sp-sc_ranked_markers.rds")
)

# * ptcellMeta enhancement
# add cl, supertype, subclass
pt_nn_tfsum <- readRDS(
  file.path(projd, "03.integration", "out",
    "tfnn_k8_cca_k5", "sum_tfnn_k8_cca_k5.seurat.rds"))
pt_neu_tfsum <- readRDS(
  file.path(projd, "03.integration", "out",
    "tfneu_k8_cca_k5", "sum_tfneu_k8_cca_k5.seurat.rds")
)

pt_nn_tfmeta <- pt_nn_tfsum@meta.data |>
  x => x[!x$isRef, ] |>
  x => `rownames<-`(x, x$barcode)
pt_neu_tfmeta <- pt_neu_tfsum@meta.data |>
  x => x[!x$isRef, ] |>
  x => `rownames<-`(x, x$barcode)
pt_tfmeta <- rbind(pt_nn_tfmeta, pt_neu_tfmeta)

L5_to_sp <- pt_tfmeta[, c("L5", "supertype_id_label")] |>
  unique() |>
  x => `rownames<-`(x, paste0("L5-", x$L5))
L5_to_sc <- pt_tfmeta[, c("L5", "subclass_id_label")] |>
  unique() |>
  x => `rownames<-`(x, paste0("L5-", x$L5))

L5_to_allen <- merge(L5_to_sp, L5_to_sc, by = "L5")
rownames(L5_to_allen) <- paste0("L5-", L5_to_allen$L5)

ptcellMeta <- file.path(
  projd, "meta",
  "pairedtag.cell.meta.all.csv"
) |>
  data.table::fread(
    file = _, sep = ",", header = TRUE,
    data.table = FALSE
  ) |>
  x => `rownames<-`(x, x$barcode)
ptcellMeta$cl <- NA
b <- intersect(ptcellMeta$barcode, pt_tfmeta$barcode)
ptcellMeta[b, "cl"] <- pt_tfmeta[b, "cl"]

for (i in c("supertype_id_label", "subclass_id_label")) {
  ptcellMeta[[i]] <- L5_to_allen[paste0("L5-", ptcellMeta$L1_2_3_4_5), i]
}
data.table::fwrite(
  x = ptcellMeta, file = file.path(projd, "meta",
    "pairedtag.cell.meta.all.with.init.tf.csv"),
  sep = ",", col.names = TRUE, row.names = FALSE
)

# * L5 to top ranked supertypes and subclasses
nn_cnss <- readRDS(
  file.path(projd, "03.integration", "out",
    "tfnn_k8_cca_k5", "sum_tfnn.consensus.mat.rds")
)
neu_cnss <- readRDS(
  file.path(projd, "03.integration", "out",
    "tfneu_k8_cca_k5", "sum_tfneu.consensus.mat.rds")
)

L5_size <- table(ptcellMeta$L1_2_3_4_5) |>
  as.data.frame(stringsAsFactors = FALSE)|>
  setNames(nm = c("L5", "size")) |>
  x => `rownames<-`(x, paste0("L5-", x$L5))

get_top_sc_groupby_sp <- function(cnss, L5_array,
                                  onAllenLevel, top) {
  cl_score_mat <- cnss[ , L5_array]
  nL5 <- L5_size[L5_array, "size"]
  total_size <- sum(nL5)
  cl_score <- if(length(L5_array) < 2) {
    round(cl_score_mat * nL5, 3)
  } else {
    cl_score_mat %*% nL5 |>
      x => round(x[ , 1], 3)
  }
  cl_score |>
    x => data.frame(
      sc = allenclMeta[names(x), onAllenLevel],
      score = x) |>
    group_by(sc) |>
    summarise(s = round(100 * sum(score) / total_size, 3)) |>
    x => x[order(x$s, decreasing = TRUE), ] |>
    x => head(x, top) |>
    x => paste(x$sc, x$s, sep = ":", collapse = ";")
}


cnss_cl_rows <- union(
  rownames(nn_cnss$cl), rownames(neu_cnss$cl)
)
cnss_cl_cols <- union(
  colnames(nn_cnss$cl), colnames(neu_cnss$cl)
)

cnss_cl <- matrix(0,
  nrow = length(cnss_cl_rows),
  ncol = length(cnss_cl_cols),
  dimnames = list(
    cnss_cl_rows,
    cnss_cl_cols
  )
)
cnss_cl[rownames(nn_cnss$cl), colnames(nn_cnss$cl)] <- nn_cnss$cl
cnss_cl[rownames(neu_cnss$cl),colnames(neu_cnss$cl)] <- neu_cnss$cl

usp <- unique(L5_to_sp$supertype_id_label)
L5sp2topsp <- vapply(usp, \(i) {
  with(L5_to_sp, paste0("L5-", L5[supertype_id_label == i])) |>
    get_top_sc_groupby_sp(cnss = cnss_cl, L5_array = _,
      onAllenLevel = "supertype_id_label",top = 3)
}, "a:30")

L5sp2topsp_df <- data.frame(
  L5sp = usp,
  top3sp = L5sp2topsp
)

L5sp2topsc <- vapply(usp, \(i) {
  with(L5_to_sp, paste0("L5-", L5[supertype_id_label == i])) |>
    get_top_sc_groupby_sp(cnss = cnss_cl, L5_array = _,
      onAllenLevel = "subclass_id_label", top = 3)
}, "a:30")
L5sp2topsc_df <- data.frame(
  L5sp = usp,
  top3sc = L5sp2topsc
)

L5sp_topsp_topsc <- merge(L5sp2topsp_df, L5sp2topsc_df, by = "L5sp")
data.table::fwrite(
  x = L5sp_topsp_topsc,
  file = file.path(projd, "meta", "L5sp_to_top_supertype_subclass.tf.init.csv"),
  sep = ",", col.names = TRUE, row.names = FALSE
)

data.table::fwrite(
  x = L5_to_sp,
  file = file.path(projd, "meta", "L5_to_sp.tf.init.csv"),
  sep = ",", col.names = TRUE, row.names = FALSE
)

# * provide cell type proportions per region
# - allen
allencellMeta <- file.path(
  projd, "meta",
  "allen10xv3.meta.with.ptregion.csv"
) |>
  data.table::fread(file = _, sep = ",", header = TRUE,
    data.table = FALSE)
allencellMeta <- subset(allencellMeta, ptr != "none")

allenclMeta <- data.table::fread(
  file = file.path(projd, "meta", "AIT21_annotation_freeze_081523.tsv"),
  sep = "\t", header = TRUE, data.table = FALSE
) |>
  x => `rownames<-`(x, paste0(clprefix, x$cl))

allencellMeta$sp <- allenclMeta[paste0(clprefix, allencellMeta$cl),
  "supertype_id_label"]
allencellMeta$sc <- allenclMeta[paste0(clprefix, allencellMeta$cl),
  "subclass_id_label"]


get_ptr2allen_celltype_mat <- function(onlevel = "sp") {
  all_sps <- allencellMeta[[onlevel]]
  all_ptr <- allencellMeta$ptr
  uptrs <- unique(all_ptr)
  asps <- unique(all_sps)

  ptr2asp_mat <- matrix(0,
    nrow = length(uptrs), ncol = length(asps),
    dimnames = list(
      uptrs,
      asps
    ))
  for (i in uptrs) {
    s <- all_sps[all_ptr == i]
    ss <- table(s)
    ptr2asp_mat[, names(ss)] <- ss
  }
  saveRDS(ptr2asp_mat,
    file = file.path(projd, "03.integration", "out",
      str_glue("allen_ptregion2{onlevel}_mat.rds")))
  return(ptr2asp_mat)
}

ptr2asp_mat <- get_ptr2allen_celltype_mat(onlevel = "sp")
ptr2asc_mat <- get_ptr2allen_celltype_mat(onlevel = "sc")

# - pairedtag
ptcellMeta <- data.table::fread(
  file = file.path(projd, "meta",
    "pairedtag.cell.meta.all.with.init.tf.csv"),
  sep = ",", header = TRUE,
  data.table = FALSE
)

all_sps <- ptcellMeta$supertype_id_label
usps <- unique(all_sps)
all_ptr <- ptcellMeta$brainregion
uptrs <- unique(all_ptr)
ptr2ptsp_mat <- matrix(0,
  nrow = length(uptrs), ncol = length(usps),
  dimnames = list(
    uptrs,
    usps
  ))

for (i in uptrs) {
  s <- all_sps[all_ptr == i]
  ss <- table(s)
  ptr2ptsp_mat[, names(ss)] <- ss
}

saveRDS(ptr2ptsp_mat,
  file = file.path(projd, "03.integration", "out",
    "pt_ptregion2sp_mat.rds"))

a <- readRDS(file.path(projd, "03.integration", "out",
  "pt_ptregion2sp_mat.rds"))

b <-  readRDS(file.path(projd, "03.integration", "out",
  "allen_ptregion2sc_mat.rds"))
