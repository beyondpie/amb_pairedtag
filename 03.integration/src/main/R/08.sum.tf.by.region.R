library(Seurat)
library(tidyverse)

# * meta
projd <- here::here()
workd <- file.path(projd, "03.integration")
mrs <- c("AMY", "CPU", "HYP", "HIP", "ERC", "NAC", "VTA", "PFC")
byf <- "vf"
tfneu_region_dir <- file.path(
  workd, "out",
  str_glue("tfneu_{byf}_region_cca_k5")
)

allenclMeta <- data.table::fread(
  file = file.path(projd, "meta", "AIT21_annotation_freeze_081523.tsv"),
  header = TRUE,
  data.table = FALSE
) |>
  x => `rownames<-`(x, paste0("cl-", x$cl))

# L5 per region mapping to super type
load_sumtf_seu <- function(r) {
  readRDS(file.path(tfneu_region_dir,
    str_glue("sumtfneu_{byf}_{r}_cca_k5.seu.rds")))
}
tfSeus <- lapply(mrs, load_sumtf_seu) |>
  setNames(object = _, nm = mrs)

# save single-cell level mapping
barcodeMetaList <- lapply(tfSeus, \(i) {
  m <- i@meta.data
  m <- m[!m$isRef, c("barcode", "cl")] |>
    as.data.frame()
  m$cl <- paste0("cl-", m$cl)
  return(m)
})

barcode2cl <- do.call(rbind, barcodeMetaList)
data.table::fwrite(x = barcode2cl,
  file = file.path(projd, "03.integration", "out",
    "tfneu_vf_region_cca_k5", "pt.raw.barcode2cl.tfbyregion.csv"),
  sep = ",", col.names = TRUE, row.names = FALSE)


get_mode_cl <- function(x) {
  names(which.max(table(x)))
}

getL5tocl <- function(seu) {
  m <- seu@meta.data
  m <- m[!m$isRef, c("barcode", "cl", "L5")] |>
    group_by(L5) |>
    summarise(cl = get_mode_cl(cl),
      ncell = length(L5)) |>
    as.data.frame() |>
    x => `rownames<-`(x, paste0("L5-", x$L5))
  m$cl <- paste0("cl-", m$cl)
  m$sp <- allenclMeta[m$cl, "supertype_label"]
  m$sc <- allenclMeta[m$cl, "subclass_label"]
  m$class <- allenclMeta[m$cl, "class_label"]
  return(m)
}

L5toCls <- lapply(tfSeus, getL5tocl) |>
  setNames(object = _, nm = mrs)

allL5s <- do.call(
  what = "c",
  args = lapply(L5toCls, \(i) {
    i$L5
  })
) |>
  unique() |>
  sort()

allL5strs <- paste0("L5-", allL5s)

## L5toCl_unsort <- lapply(mrs, \(region) {
##   r <- rep("NA", length(allL5strs)) |>
##     setNames(object = _, nm = allL5strs)
##   m <- L5toCls[[region]]
##   r[rownames(m)] <- str_glue("{region}:{m$sp}@{m$ncell}")
##   return(r)
## }) |>
##   x => do.call(what = "rbind", args = x)

## L5toCl_sort <- vapply(seq_along(allL5strs), \(i) {
##   raw <- L5toCl_unsort[, i]
##   ns <- vapply(L5toCl_unsort[, i], \(j) {
##     as.integer(str_split_1(j, pattern = "@")[2])
##   }, 1L)
##   paste(raw[order(ns, decreasing = TRUE)], collapse = ";") |>
##     gsub("NA", "", x = _) |>
##     gsub(";+", ";", x = _) |>
##     gsub(";$", "", x = _)
## }, "c")

## L5toCl_df <- data.frame(
##   L5 = allL5s,
##   tfbyregion = L5toCl_sort
## )

## data.table::fwrite(
##   x = L5toCl_df,
##   file = file.path(tfneu_region_dir,
##     "L5tosupertype_tfregion_sortbyncell.csv"),
##   sep = ",",
##   col.names = TRUE, row.names = FALSE
## )

# * check cl mode
## r <- "AMY"
## a <- readRDS(file.path(
##   tfneu_region_dir,
##   str_glue("tf_{r}/query.with.tf-cca-kac5_on-cl.rds")
## ))
## # here cl is string
## b2pcl <- data.frame(
##   barcode = a$barcode,
##   L5 = a$L1_2_3_4_5,
##   cl = a$predicted.id
## )
# cl is factor in this data
## a <- readRDS(file.path(projd, "data",
##   "allen_seurat", "allen.10xv3.AMY.neu.clean.all.rds"))
## cls <- a$cl
## # pasting factor will use original string that factor uses.
## cl_cls <- paste0("cl_", cls)
## b <- factor(c(rep("a", 1000), rep("b", 100)))
## str(paste0("factor-", b))

# * add global integration result and the real cell numbers per region
neu_ptseu <- readRDS(
  file.path(workd, "out", "tfneu_k8_cca_k5",
    "sum_tfneu_k8_cca_k5.seurat.rds"))
L5tosp <- neu_ptseu@meta.data |>
  x => x[!x$isRef, ] |>
  group_by(L5) |>
  summarise(sp = unique(supertype_id_label)) |>
  as.data.frame() |>
  x => x[, c("L5", "sp")] |>
  x => `rownames<-`(x, paste0("L5-", x$L5))

ptMeta <- data.table::fread(
  file = file.path(projd, "meta", "pairedtag.cell.meta.all.csv"),
  header = TRUE,
  sep = ",",
  data.table = FALSE
)

pt2mr <- data.frame(
  pt = c("HYP", "CPU", "HCa", "HCp", "ERC", "AMY", "NAC", "VTA_SnR", "PFC"),
  mr = c("HYP", "CPU", "HIP", "HIP", "ERC", "AMY", "NAC", "VTA", "PFC")
) |>
  x => `rownames<-`(x, x$pt)

get_L5_mr_count <- function(vec_of_region) {
  t <- pt2mr[vec_of_region, "mr"] |>
    table()
  paste(names(t), t, sep = ":", collapse = ";")
}

nmr_to_vec <- function(str_of_nmr) {
  r <- rep("0", length(mrs)) |>
    setNames(object = _, nm = mrs)
  a <- str_split_1(str_of_nmr, ";")
  b <- vapply(a, str_split_1, c("AMY", "282"), pattern = ":")
  bb <- b[2, ] |>
    setNames(object = _, nm = b[1, ])
  r[names(bb)] <- bb
  return(r)
}

L5tomr_cnt <- ptMeta[ptMeta$L1_2_3_4_5 %in% allL5s, ] |>
  group_by(L1_2_3_4_5) |>
  summarise(
    mcnt = get_L5_mr_count(brainregion)) |>
  as.data.frame() |>
  x => `rownames<-`(x, paste0("L5-", x$L1_2_3_4_5))

L5tomr_cnt_list <- lapply(rownames(L5tomr_cnt), \(l5){
  nmr_to_vec(L5tomr_cnt[l5, "mcnt"])
}) |>
  setNames(object = _, nm = rownames(L5tomr_cnt))

# * update L5tomr_cnt
allL5strs <- paste0("L5-", allL5s)
L5toCl_unsort <- lapply(mrs, \(region) {
  r <- rep("NULL", length(allL5strs)) |>
    setNames(object = _, nm = allL5strs)
  m <- L5toCls[[region]]
  ntotal <- vapply(rownames(m), \(l5) {
    L5tomr_cnt_list[[l5]][region]
  }, "20") 
  r[rownames(m)] <- str_glue("{region}:{m$sp}@{m$ncell}@{ntotal}")
  return(r)
}) |>
  x => do.call(what = "rbind", args = x)

L5toCl_sort <- vapply(seq_along(allL5strs), \(i) {
  raw <- L5toCl_unsort[, i]
  ns <- vapply(L5toCl_unsort[, i], \(j) {
    as.integer(str_split_1(j, pattern = "@")[2])
  }, 1L)
  paste(raw[order(ns, decreasing = TRUE)], collapse = ";") |>
    gsub("NULL", "", x = _) |>
    gsub(";+", ";", x = _) |>
    gsub(";$", "", x = _)
}, "c")

tfbyglobal <- gsub("^\\d+ ", "",L5tosp[allL5strs, "sp"])
L5toCl_df <- data.frame(
  L5 = allL5s,
  tfbyregion = L5toCl_sort,
  tfbyglobal = tfbyglobal
)

top1_sc_by_region <- vapply(L5toCl_df$tfbyregion, \(l) {
  ll <- str_split_1(l, ";")[1] |>
    x => str_split_1(string = x, pattern = "@")[1] |>
    x => str_split_1(string = x, pattern = ":")[2] |>
    gsub("_\\d", "", x = _)
}, "str")
sc_by_global <- gsub("_\\d", "", x = tfbyglobal)

a <- sc_by_global == top1_sc_by_region
aa<- rep("FALSE", length(allL5s))
names(aa) <- allL5strs
aa[a] <- "TRUE"
L5toCl_df <- data.frame(
  L5 = allL5s,
  tfbyregion = L5toCl_sort,
  tfbyglobal = tfbyglobal,
  matchby_subclass = aa
)
# about 2/3 of them are aligned on subclass-level
# how about the unaligned one: their size

data.table::fwrite(
  x = L5toCl_df,
  file = file.path(tfneu_region_dir,
    "L5tosp_byregion_byglobal_match.csv"),
  sep = ",",
  col.names = TRUE, row.names = FALSE
)

# * hard threshold
L5tomr_cnt_df <- do.call(rbind, L5tomr_cnt_list) |>
  apply(X = _, MARGIN = 2, as.numeric) |>
  x => `rownames<-`(x, names(L5tomr_cnt_list))

get_total_ncell_after_remove <- function(n) {
  r <- L5tomr_cnt_df
  r[r <= n] <- 0
  sum(r)
}

total <- sum(L5tomr_cnt_df)
get_total_ncell_after_remove(n = 10) * 100/ total
get_total_ncell_after_remove(n = 20) * 100/ total
total-get_total_ncell_after_remove(n = 10 ) 
get_total_ncell_after_remove(n = 30) * 100/ total

x <- seq(from = 1, to = 9) * 10
y <- vapply(x, \(i){
  round(get_total_ncell_after_remove(n = i) * 100 / total,3)},
  0.99)
L5_size <- rowSums(L5tomr_cnt_df)
L5toCl_df$ntotal <- L5_size[rownames(L5toCl_df)]
with(L5toCl_df, quantile(ntotal[matchby_subclass == "FALSE"]))
with(L5toCl_df, quantile(ntotal[matchby_subclass == "TRUE"]))


# * use 10 as threshold per region to filter cells per L5.
hard_threshold <- 10
allL5strs <- paste0("L5-", allL5s)
L5toCl_filter <- lapply(mrs, \(region) {
  r <- rep("NULL", length(allL5strs)) |>
    setNames(object = _, nm = allL5strs)
  m <- L5toCls[[region]]
  ntotal <- vapply(rownames(m), \(l5){
    as.integer(L5tomr_cnt_list[[l5]][region])
  }, 20L)
  m <- m[ntotal > hard_threshold, ]
  r[rownames(m)] <- m$sp
  return(r)}) |>
  x => do.call(what = "rbind", args = x) |>
  x => `rownames<-`(x, mrs)
L5toCl_filter_tri <- reshape2::melt(L5toCl_filter) |>
  setNames(object = _, nm = c("region", "L5str", "sp")) |>
  subset(x = _, subset = sp != "NULL")
L5toCl_filter_tri$region <- levels(
  L5toCl_filter_tri$region)[L5toCl_filter_tri$region]
L5toCl_filter_tri$L5str <- levels(
  L5toCl_filter_tri$L5str)[L5toCl_filter_tri$L5str]
L5toCl_filter_tri$L5r <- with(
  L5toCl_filter_tri, paste(L5str, region, sep = ":")
)
rownames(L5toCl_filter_tri) <- L5toCl_filter_tri$L5r

# barcode to L5r
ptMeta <- data.table::fread(
  file = file.path(projd, "meta", "pairedtag.cell.meta.all.csv"),
  header = TRUE,
  sep = ",",
  data.table = FALSE
)
ptMeta$majoregion <- pt2mr[ptMeta$brainregion, "mr"]
ptMeta$L5r <- paste0("L5-", ptMeta$L1_2_3_4_5)
rownames(ptMeta) <- ptMeta$barcode

psL5r <- with(ptMeta, paste(L5r, majoregion, sep = ":"))
index <- psL5r %in% L5toCl_filter_tri$L5r

ptMeta$L5rsp <- "NULL"
ptMeta$L5rsp[index] <-
  L5toCl_filter_tri[psL5r[index], "sp"]

ptMeta$isNeuFromL1 <- FALSE
ptMeta$isNeuFromL1[
  paste0("L5-", ptMeta$L1_2_3_4_5) %in% colnames(L5toCl_filter)
] <- TRUE

ptMeta$L5r[index] <- psL5r[index]
ptMeta$L5r[ptMeta$isNeuFromL1 & (!index)] <- "LQ_tf_by_region"
data.table::fwrite(
  x = ptMeta, file = file.path(workd, "out", "tfneu_vf_region_cca_k5",
    "pairedtag.cell.meta.all.with.neutfbyregion.v3.csv"),
  sep = ",",
  row.names = FALSE, col.names = TRUE
)

# * [Optional] plot transfer label score
# * add non-neuronal cell labels
outfnm <- file.path(projd, "meta",
  "pairedtag.cell.meta.all.with.tfv3.csv")

ptMeta_neuLabel <- data.table::fread(
  file = file.path(workd, "out", "tfneu_vf_region_cca_k5",
    "pairedtag.cell.meta.all.with.neutfbyregion.v3.csv"),
  sep = ",",
  header = TRUE, data.table = FALSE
) |>
  x => `rownames<-`(x, x$barcode)

ptMeta_inittf <- data.table::fread(
  file = file.path(projd, "meta",
    "pairedtag.cell.meta.all.with.init.tf.csv"),
  header = TRUE,
  sep = ",", data.table = FALSE
) |> x => `rownames<-`(x, x$barcode)

ptMeta <- ptMeta_inittf
ptMeta$L5r <- ptMeta_neuLabel[rownames(ptMeta), "L5r"]
ptMeta$isNeuFromL1 <- ptMeta_neuLabel[rownames(ptMeta), "isNeuFromL1"]

# now add the cls for previous region-specific mapping
L5rtoCl_byregion <- lapply(mrs, \(mr) {
  a <- L5toCls[[mr]]
  data.frame(
    L5r = paste0("L5-", a$L5, ":", mr),
    cl = a$cl
  )
}) |> do.call(rbind, args = _) |>
  x => `rownames<-`(x, x$L5r)

index <- ptMeta$L5r %in% rownames(L5rtoCl_byregion)
ptMeta$cl[index] <- L5rtoCl_byregion[ptMeta$L5r[index], "cl"]

nn_L5rtocl <- ptMeta[!ptMeta$isNeuFromL1, ] |>
  x => x[!is.na(x$cl), ] |>
  group_by(L5r) |>
  summarise(ucl = get_mode_cl(cl)) |>
  as.data.frame() |>
  x => `rownames<-`(x, x$L5r)

index <- ptMeta$L5r %in% rownames(nn_L5rtocl)
ptMeta$cl[index] <- paste0("cl-", nn_L5rtocl[ptMeta$L5r[index], "ucl"])

# now add L5rsp, L5rsc, L5rclass
index <- ptMeta$cl %in% rownames(allenclMeta)
for (i in c(
  "supertype_label",
  "supertype_id_label", "subclass_id_label", "class_id_label")) {
  ptMeta[[i]] <- NA
  ptMeta[[i]][index] <- allenclMeta[ptMeta$cl[index], i]
}

## check results
neuL5r <- with(ptMeta, L5r[grepl(":", L5r)])
a <- ptMeta[ptMeta$L5r %in% neuL5r, ]
all(a$subclass_label == ptMeta_neuLabel[rownames(a), "L5rsp"])
nnL5r <- setdiff(ptMeta$L5r, neuL5r) |>
  x => x[!grepl("LQ_", x)]
a <- ptMeta[ptMeta$L5r %in% nnL5r,]
all(a$supertype_id_label == ptMeta_inittf[rownames(a), "supertype_id_label"])

# remove this since it's used for checking only
ptMeta$supertype_label <- NULL

data.table::fwrite(
  x = ptMeta, file = outfnm, row.names = FALSE, col.names = TRUE
)
