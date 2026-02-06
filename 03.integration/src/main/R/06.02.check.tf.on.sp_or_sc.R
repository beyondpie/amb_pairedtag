library(tidyverse)
library(Seurat)

# shared meta
projd <- here::here()
clprefix <- "cl-"
onAllenLevel <- "supertype_id_label"

# parameters
# * neuron
tfSumdir <- file.path(projd, "03.integration", "out",
  "tfneu_k8_cca_k5")
tfsumfnm <- file.path(tfSumdir, "sum_tfneu_k8_cca_k5.seurat.rds")
cnssfnm <- file.path(tfSumdir, "sum_tfneu.consensus.mat.rds")
ptseu <- readRDS(file.path(projd, "data/pairedtag_seurat",
  "ptRNA.neu.k8.L5ds50.rds"))
onAllenLevel <- "subclass_id_label"
outfnm <- file.path(tfSumdir, str_glue("check_tfneu_k8_cca_k5.{onAllenLevel}.csv"))

# * non-neurons
## tfSumdir <- file.path(projd, "03.integration", "out",
##   "tfnn_k8_cca_k5")
## tfsumfnm <- file.path(tfSumdir, "sum_tfnn_k8_cca_k5.seurat.rds")
## cnssfnm <- file.path(tfSumdir, "sum_tfnn.consensus.mat.rds")
## ptseu <- readRDS(file.path(projd, "data/pairedtag_seurat",
##   "ptRNA.nn.k8.L5ds30.rds"))
## outfnm <- file.path(tfSumdir,
##   str_glue("check.tfnn_k8_cca_k5.{onAllenLevel}.csv"))


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

# 1. meta data needed
common_markers <- data.table::fread(
  file = file.path(projd, "meta", "common_gene_markers.csv"),
  sep = ",", header = TRUE, data.table = FALSE
) |>
  x => `rownames<-`(x, x$gene)

# pairted-tag pseudobulk-level seurat
# on supertype level
pt_pseu_seu <- readRDS(file.path(
  projd, "03.integration", "out",
  "pt.pseudo_seurat.supertype.rds"
))

allenclMetafnm <- file.path(projd, "meta",
  "AIT21_annotation_freeze_081523.tsv")
allenclMeta <- data.table::fread(
  file = allenclMetafnm, sep = "\t", header = TRUE,
  data.table = FALSE
)
rownames(allenclMeta) <- paste0(clprefix, allenclMeta$cl)
# filter by subclass we focus on
allenclMeta <- allenclMeta[allenclMeta[[onAllenLevel]] %in% tfscs, ]

allencellMetafnm <- file.path(projd, "data/allen",
  "allen.10xv3.cell.meta.csv")
allencellMeta <- data.table::fread(
  file = allencellMetafnm, sep = ",", header = TRUE,
  data.table = FALSE)
rownames(allencellMeta) <- allencellMeta$barcode
# filter by subclass we focus on
allencellMeta <- subset(allencellMeta, allencellMeta$cl %in% allenclMeta$cl)
allencellMeta[[onAllenLevel]] <- allenclMeta[
  paste0(clprefix, allencellMeta$cl), onAllenLevel]

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

# 2. allen region to ours
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

a2pt <- lapply(aregions, \(a) {
  with(pt2aRegion, pt[allen %in% a])
}) |> setNames(object = _, nm = aregions)


# 3. get allen subclass meta data
# - subclass gene marker
rank_markers <- function(mgs) {
  g <- str_split_1(mgs, ",")
  table(g) |>
    as.data.frame(x = _, stringsAsFactors = FALSE) |>
    x => `rownames<-`(x, x$genes) |>
    x => x[order(x$Freq, decreasing = TRUE), ]
}

get_topk_markers <- function(mgs, top = 5) {
  rank_markers(mgs) |>
    head(x = _, n = top)
}


sc2mgs <- allenclMeta |>
  group_by(.data[[onAllenLevel]]) |>
  summarise(marker = paste(cluster.markers, collapse = ",")) |>
  as.data.frame()

sc2top_markers <- lapply(seq_len(nrow(sc2mgs)), function(i) {
  get_topk_markers(sc2mgs[i, "marker"])
}) |>
  setNames(object = _, nm = sc2mgs[[onAllenLevel]])

# - subclass major regions (paired-tag regions)
rank_regions <- function(roi, top = 3) {
  gsub(" |  ", "", x = roi) |>
    gsub("-", "_", x = _) |>
    table() |>
    as.data.frame(x = _, stringsAsFactors = FALSE) |>
    x => x[order(x$Freq, decreasing = TRUE), ] |>
    mutate(ratio = round(Freq / sum(Freq), digits = 3) * 100)
}

get_top_regions <- function(roi_stat, top = 3) {
  roi_stat |>
     x => head(x, n = top) |>
    x => with(x, paste(Var1, ratio, sep = ":", collapse = ";"))
}

get_top_pt_region <- function(roi_stat, top = 3) {
  n_total <- sum(roi_stat$Freq)
  r_df <- vapply(seq_len(nrow(roi_stat)), \(i) {
    r <- rep(0.0, length(ptregions)) |>
      setNames(object = _, nm = ptregions)
    a <- roi_stat[i, 1]
    if (a %in% names(a2pt)) {
      r[a2pt[[a]]] <- roi_stat[i, 2]
    }
    return(r)
  }, rep(0.0, length(ptregions)))
  r <- round(rowSums(r_df) / n_total, 3) * 100
  r |>
    x => x[order(x, decreasing = TRUE)] |> 
    x => x[x >= 0.1] |> 
    x => head(x, n = top) |>
    x => paste(names(x), x, sep = ":", collapse = ";")
}

## roi <- allencellMeta$roi[
##   allencellMeta$subclass_id_label %in% allenclMeta$subclass_id_label[1]]
## roi_stat <- rank_regions(roi)

sc2rgs <- allencellMeta |>
  group_by(.data[[onAllenLevel]]) |>
  summarise(
    allen_region = get_top_regions(rank_regions(roi), top = 3),
    pt_region = get_top_pt_region(rank_regions(roi), top = 3)
  ) |>
  as.data.frame() |>
  subset(x = _, nchar(pt_region) > 1) |>
  x => `rownames<-`(x, x[[onAllenLevel]])

# - sex ratio
sc2sex <- allenclMeta |>
  group_by(.data[[onAllenLevel]]) |>
  summarise(male = round(mean(M), 3) * 100) |>
  as.data.frame() |>
  x => `rownames<-`(x, x[[onAllenLevel]])

# - 10xv3 ratio

sc2nv3 <- allenclMeta |>
  group_by(.data[[onAllenLevel]]) |>
  summarise(nv3 = sum(v3.size)) |>
  as.data.frame() |>
  x => `rownames<-`(x, x[[onAllenLevel]])

sc2v3 <- allenclMeta |>
  group_by(.data[[onAllenLevel]]) |>
  summarise(v3 =
              round(sum(v3.size) * 100 /
                      (sum(v3.size) + sum(v2.size) + sum(multiome.size)), 3)) |>
  as.data.frame() |>
  x => `rownames<-`(x, x[[onAllenLevel]])

# 4. group our L5-cluster to ptsp based on super-type
L5_sc <- tfsum@meta.data[!tfsum$isRef, ] |>
  group_by(L5) |>
  summarise(sc = unique(.data[[onAllenLevel]])) |>
  as.data.frame() |>
  x => `rownames<-`(x, paste0("L5-", x$L5))
L5_sp <- tfsum@meta.data[!tfsum$isRef, ] |>
  group_by(L5) |>
  summarise(sp = unique(supertype_id_label)) |>
  as.data.frame() |>
  x => `rownames<-`(x, paste0("L5-", x$L5))

get_mode_cl <- function(x) {
  names(which.max(table(x)))
}

L5_cl <- tfsum@meta.data[!tfsum$isRef, ] |>
  group_by(L5) |>
  summarise(cl = get_mode_cl(cl)) |>
  as.data.frame() |>
  x => `rownames<-`(x, paste0("L5-", x$L5))

# 5. get ptsp meta data
# - get brain region
ptcellMeta$supertype_id_label <- L5_sp[
  paste0("L5-", ptcellMeta$L1_2_3_4_5), "sp"]
ptcellMeta$cl <- L5_cl[
  paste0("L5-", ptcellMeta$L1_2_3_4_5), "cl"]

all_sp <- unique(ptcellMeta$supertype_id_label)

L5sp_rgs <- ptcellMeta |>
  group_by(supertype_id_label) |>
  summarize(pt_region = get_top_regions(rank_regions(brainregion),
    top = 3
  )) |>
  as.data.frame() |>
  x => `rownames<-`(x, x$supertype_id_label)

# - get sex raio
L5sp_male <- ptcellMeta |>
  group_by(supertype_id_label) |>
  summarize(pt_male = round(sum(sex == "Male") / length(sex), 3) * 100) |>
  as.data.frame() |>
  x => `rownames<-`(x, x$supertype_id_label)
  
# - get replicate ratio
get_pt_rep <- function(rep_array) {
  all_repA <- 100 * round(
    sum(grepl("A", rep_array)) / length(rep_array), 3)
  male_repA <- 100 * round(
    sum(grepl("MaleA", rep_array)) / sum(grepl("Male", rep_array)), 3)
  female_repA <- 100 * round(
    sum(grepl("FemaleA", rep_array)) / sum(grepl("Female", rep_array)), 3)
  str_glue("allA:{all_repA};mA:{male_repA};fA:{female_repA}")
}

L5sp_rep <- ptcellMeta |>
  group_by(supertype_id_label) |>
  summarize(pt_rep = get_pt_rep(rep)) |>
  as.data.frame() |>
  x => `rownames<-`(x, x$supertype_id_label)

# - get histone ratio
L5sp_mod <- ptcellMeta |>
  group_by(supertype_id_label) |>
  summarize(pt_mod = get_top_regions(rank_regions(modality),top = 4)) |>
  as.data.frame() |>
  x => `rownames<-`(x, x$supertype_id_label)

# 6. get top ranked subclasses the ptsp mapped to.
L5sp_nL5 <- ptcellMeta |>
  group_by(supertype_id_label) |>
  summarize(nL5 = length(unique(L1_2_3_4_5))) |>
  as.data.frame() |>
  x => `rownames<-`(x, x$supertype_id_label)

L5_size <- table(ptcellMeta$L1_2_3_4_5) |>
  as.data.frame(stringsAsFactors = FALSE)|>
  setNames(nm = c("L5", "size")) |>
  x => `rownames<-`(x, paste0("L5-", x$L5))

L5sp_size <- table(ptcellMeta$supertype_id_label) |>
  as.data.frame(stringsAsFactors = FALSE) |>
  setNames(nm = c("sp", "ncell")) |>
  x => `rownames<-`(x, x$sp)

# for test the function below
# L5_array <- L5_sp$L5[L5_sp$sp == "0008 L5/6 IT TPE-ENT Glut_2"]
get_top_sc_groupby_sp <- function(L5_array, top = 3 ) {
  cl_score_mat <- cnss$cl[ , paste0("L5-", L5_array)]
  nL5 <- L5_size[paste0("L5-", L5_array), "size"]
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
    x => head(x, top)
}

print_top_sc <- function(L5_array, top = 3) {
  get_top_sc_groupby_sp(L5_array, top = top) |>
    x => paste(x$sc, x$s, sep = ":", collapse = ";")
}

L5sp_top_scs <- L5_sp |>
  group_by(sp) |>
  summarize(top_scs = print_top_sc(L5, top = 3))

L5sp2topscs <- lapply(all_sp, \(i) {
  with(L5_sp, L5[sp == i]) |>
    get_top_sc_groupby_sp(L5_array = _, top = 3)
})
names(L5sp2topscs) <- all_sp
L5sp_ntopscs <- vapply(L5sp2topscs, nrow, 1)

# 7. link the info to one Excel.
sp_basic_info <- data.frame(
  name = all_sp,
  ncell = L5sp_size[all_sp, 2],
  ncl = L5sp_nL5[all_sp, 2],
  rep = L5sp_rep[all_sp, 2],
  male = L5sp_male[all_sp, 2],
  mod = L5sp_mod[all_sp, 2],
  region = L5sp_rgs[all_sp, 2]
) |>
  x => `rownames<-`(x, x$name)

get_topk_sc <- function(sps, k = 1) {
  vapply(sps, \(i) {
    L5sp2topscs[[i]][k, "sc", drop = TRUE]
  }, "001 allen subclass")
}
get_topk_sc_score <- function(sps, k = 1) {
  vapply(sps, \(i) {
    L5sp2topscs[[i]][k, "s", drop = TRUE]
  }, 100.0)
}

get_sc_info <- function(scs) {
  data.frame(
    ncellv3 = sc2nv3[scs, 2],
    nv3ratio = sc2v3[scs, 2],
    male = sc2sex[scs, 2],
    allen_region = sc2rgs[scs, "allen_region"],
    pt_region = sc2rgs[scs, "pt_region"]
  )
}

top1_sc <- get_topk_sc(all_sp, k = 1)
top2_sc <- get_topk_sc(all_sp, k = 2)
top3_sc <- get_topk_sc(all_sp, k =3)
top1_sc_info <- top1_sc |>
  get_sc_info() |>
  x => `rownames<-`(x, all_sp)
top2_sc_info <- top2_sc |>
  get_sc_info() |>
  x => `rownames<-`(x, all_sp)
top3_sc_info <- top3_sc |>
  get_sc_info() |>
  x => `rownames<-`(x, all_sp)

# get gene expression
tf_ptMeta <- tfsum@meta.data[!tfsum$isRef, ] |>
  x => `rownames<-`(x, x$barcode)

norm_exp <- pt_pseu_seu@assays$RNA$data |>
  setNames(object = _, nm = colnames(pt_pseu_seu)) |>
  x => `rownames<-`(x, rownames(pt_pseu_seu))
scale_exp <- pt_pseu_seu@assays$RNA$scale.data |>
  setNames(object = _, nm = colnames(pt_pseu_seu)) |>
  x => `rownames<-`(x, rownames(pt_pseu_seu))

get_common_exp <- function(exp_mat, sps) {
  cm <- common_markers[
    common_markers$gene %in% rownames(exp_mat), ]
  vapply(sps, \(i) {
    exps <- round(exp_mat[cm$gene,
      gsub("_", "-", paste0("g", i))], 2)
    paste(cm$label,
      paste(cm$gene, exps, sep = ":"),
      sep = "@", collapse = ";")
  }, "neu@Snap25:0.1;inh@Slcblabla")
}

sp_com_norm_exp <- data.frame(
  name = all_sp,
  classic_markers_norm_exp = get_common_exp(norm_exp, all_sp)
)
sp_com_scale_exp <- data.frame(
  name = all_sp,
  classic_markers_scale_exp = get_common_exp(scale_exp, all_sp)
)
sp_com_exp <- merge(
  sp_com_norm_exp,
  sp_com_scale_exp,
  by = "name")

sp_basic_info <- merge(sp_basic_info, sp_com_exp,by = "name")

sc2top_markers_flt <- lapply(sc2top_markers, \(i) {
  i$g[i$g %in% tfgenes]
})

sp2sc_info <- data.frame(
  name = all_sp,
  top1_sc = top1_sc,
  top1_sc_info = top1_sc_info,
  top2_sc = top2_sc,
  top2_sc_info = top2_sc_info,
  top3_sc = top3_sc,
  top3_sc_info = top3_sc_info
) |>
  x => `rownames<-`(x, all_sp)

get_marker_exp <- function(exp_mat, sps, k = 1) {
  vapply(sps, \(i) {
    sc <- sp2sc_info[i, str_glue("top{k}_sc")]
    genes <- sc2top_markers_flt[[sc]]
    exps <- round(exp_mat[genes, gsub("_", "-", paste0("g", i))], 2)
    paste(genes, exps, sep = ":", collapse = ";")
  }, "a:0.1;b:0.2")
}

top1_markers_norm_exp <- get_marker_exp(norm_exp, all_sp, k = 1)
top2_markers_norm_exp <- get_marker_exp(norm_exp, all_sp, k = 2)
top3_markers_norm_exp <- get_marker_exp(norm_exp, all_sp, k = 3)

top1_markers_scale_exp <- get_marker_exp(scale_exp, all_sp, k = 1)
top2_markers_scale_exp <- get_marker_exp(scale_exp, all_sp, k = 2)
top3_markers_scale_exp <- get_marker_exp(scale_exp, all_sp, k = 3)

sp2sc_info <- data.frame(
  name = all_sp,
  top1_sc = top1_sc,
  top1_sc_tfscore = get_topk_sc_score(all_sp, k = 1),
  top1_sc_info = top1_sc_info,
  top1_markers_norm_exp = top1_markers_norm_exp,
  top1_markers_scale_exp = top1_markers_scale_exp,
  top2_sc = top2_sc,
  top2_sc_tfscore = get_topk_sc_score(all_sp, k =2),
  top2_sc_info = top2_sc_info,
  top2_markers_norm_exp = top2_markers_norm_exp,
  top2_markers_scale_exp = top2_markers_scale_exp,
  top3_sc = top3_sc,
  top3_sc_tfscore = get_topk_sc_score(all_sp, k =3),
  top3_sc_info = top3_sc_info,
  top3_markers_norm_exp = top3_markers_norm_exp,
  top3_markers_scale_exp = top3_markers_scale_exp
) |>
  x => `rownames<-`(x, all_sp)

# merge info
sp_info <- merge(
  sp_basic_info,
  sp2sc_info,
  by = "name"
)

data.table::fwrite(
  x = sp_info,
  file = outfnm,
  sep = ",", row.names = F, col.names = T
)
