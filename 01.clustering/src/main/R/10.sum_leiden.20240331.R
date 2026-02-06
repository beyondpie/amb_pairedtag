library(tidyverse)
library(rlang)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
# inhibit scientific notation
# for cluster label.
options(scipen=999)

projroot <- here::here()
workdir <- file.path(projroot, "01.clustering")
rscdir <- file.path(workdir, "src/main/resource")

# * L1 level summary
icll <- 1
pcll <- paste0("L", icll-1)
cll <- paste0("L", icll)
rdir <- file.path(workdir, "out", "tscc2", cll)
prefix <- paste0("RNA_k8_npc50_k50_", pcll)

c_j <- 0
leidenfnm <- file.path(rdir,
  paste0(prefix, "_c", c_j, ".cellmeta.leiden.csv"))
# avgSil: 0.24, cl = 17
bestld <- "r0.3_seed3"

ldsum <- data.table::fread(file = leidenfnm, sep = ",",
  header = TRUE, data.table = FALSE)

cellmeta <- data.frame(
  barcode = ldsum$barcode,
  umap1 = ldsum$umap1,
  umap2 = ldsum$umap2
)
cellmeta[[cll]] <- ldsum[[bestld]]

data.table::fwrite(
  x = cellmeta,
  file = file.path(rscdir,
    paste0("pt.barcode2group.", cll, ".20240331.csv")),
  row.names = FALSE, col.names = TRUE
)

# * L2 level summary
icll <- 2
pcll <- paste0("L", icll-1)
cll <- paste0("L", icll)
rdir <- file.path(workdir, "out", "tscc2", cll)
prefix <- paste0("RNA_vg3000_npc30_k50_", pcll)

bestlds <- c(
  "r0.1_seed4", # 6cls
  "r0.4_seed0", # 8cls orig: r0.1_seed1 2cls
  "r0.3_seed2", # 11cls orig: r0.1_seed7 2cls
  "r0.7_seed4", # 21cls
  "r0.4_seed4", # 18cls
  "r0.2_seed1", # 5cls
  "r0.2_seed5", # 3cls orig: r0.1_seed0 2cls
  "r0.4_seed1", # 14cls orig: r0.3_seed1 11cls
  "r0.3_seed4", # 8cls orig: r0.1_seed8 5cls
  "r0.5_seed9", # 11cls
  "r0.1_seed3", # 4cls
  "r0.3_seed9", # 9cls
  "r0.1_seed3", # 4cls
  "r0.6_seed8", # 21cls orig: r0.1_seed4 2cls
  "r0.7_seed3", # 23cls
  "r0.4_seed5", # 14cls
  "r0.1_seed0" # 2cls
  )
names(bestlds) <- paste0("c", 0:16)
# read cell meta based on bestlds
cellmetaL2 <- lapply(0:16, \(c_j) {
  nm <- paste0("c", c_j)
  t <- data.table::fread(
    file.path(
      rdir,
      paste0(prefix, "_c", c_j, ".cellmeta.leiden.csv")
    ),
    sep = ",", header = TRUE, data.table = FALSE
  )
  r <- data.frame(
    barcode = t$barcode,
    L2UMAP1 = t$umap1,
    L2UMAP2 = t$umap2
  )
  r[[cll]] <- t[[bestlds[nm]]]
  rownames(r) <- r$barcode
  message(c_j, " has ", length(table(r[[cll]])), " cls.")
  return(r)
}) |> do.call("rbind", args = _)

# append L1 result
cellmetaL1 <- data.table::fread(
  file.path(rscdir, "pt.barcode2group.L1.20240331.csv"),
  sep = ",", header = TRUE, data.table = FALSE
) |>
  dplyr::rename(L1UMAP1 = umap1, L1UMAP2 = umap2) |>
  x => `rownames<-`(x, x$barcode)

cellmeta <- merge(cellmetaL1, cellmetaL2, by = "barcode")
cellmeta$L1_2 <- cellmeta$L1 * 100 + cellmeta$L2
data.table::fwrite(
  x = cellmeta,
  file = file.path(rscdir, "pt.barcode2group.L2.20240331.csv"),
  sep = ",", row.names = FALSE, col.names = TRUE
)

# * summarize L3-level clustering
# get bestlds automatically then check manually for each of them.
cellmetaL2 <- data.table::fread(
  file = file.path(rscdir, "pt.barcode2group.L2.20240331.csv"),
  sep = ",", header = TRUE, data.table = FALSE
) |> x => `rownames<-`(x, x$barcode)

icll <- 3
pcll <- paste0("L", icll-1)
cll <- paste0("L", icll)
rdir <- file.path(workdir, "out", "tscc2", cll)
prefix <- "RNA_vg2000_npc30_k30_L1_2"

getBestlds <- function(c_j) {
  getrcols <- function(r, colnms) {
    r <- format(as.numeric(r), nsmall = 1)
    ptn <- str_glue("r{format(r)}_")
    colnms[grep(ptn, colnms)]
  }
  getSeeds <- function(colnms) {
    colnms |>
      x => x[grep("seed", x)] |>
      map_chr(\(x) str_split_1(x, "_")[2]) |>
      gsub("seed", "", x = _) |>
      as.integer()
  }
  r <- seq(0.1, 2, 0.1)
  res_silht <- data.table::fread(
    file = file.path(rdir, paste0(prefix, "_c", c_j, ".silht.csv")),
    sep = ",", header = TRUE, data.table = FALSE
  )
  avgSils <- map_vec(r, \(rr) {
    mean(rowMeans(res_silht[, getrcols(rr, colnames(res_silht))]))
  }) |>
    setNames(object = _, paste0("r", format(r, nsmall = 1)))
  rr <- r[which.max(avgSils)]
  rrcols <- getrcols(rr, colnames(res_silht))
  seeds <- getSeeds(rrcols)
  ss <- seeds[which.max(colMeans(res_silht[ , rrcols]))]
  bestld <- str_glue("r{format(as.numeric(rr), nsmall=1)}_seed{ss}")
  message(str_glue(
    "c{c_j} bestld: {bestld}",
    " with sil {mean(res_silht[ , bestld])}."
  ))
  return(bestld)
}

c_js <- unique(cellmetaL2$L1_2)
bestlds <- map_chr(c_js, getBestlds) |>
  x => data.frame(
    L1_2 = paste0("c", c_js),
    bestld = x,
    note="auto"
  ) |>
  x => x[order(x$L1_2), ]

data.table::fwrite(
  x = bestlds,
  file = file.path(rscdir, "auto.bestlds.L1_2.20240401.csv"),
  row.names = FALSE, col.names = TRUE
)

# after manual check, load bestlds
bestlds <- data.table::fread(
  file = file.path(rscdir, "auto.bestlds.L1_2.20240401.csv"),
  sep = ",", header = TRUE, data.table = FALSE) |>
  x => `rownames<-`(x, x$L1_2)

cellmetaL3 <- lapply(rownames(bestlds), function(c_j) {
  fnm <- file.path(rdir, paste0(prefix, "_", c_j, ".cellmeta.leiden.csv"))
  t <- data.table::fread(fnm, sep = ",", header = TRUE, data.table = FALSE)
  bestld <- bestlds[c_j, "bestld"]
  if (! (bestld %in% colnames(t)) ) {
    stop(paste0(bestld, " not exist in ", c_j))
  }
  r <- data.frame(
    barcode = t$barcode,
    L3UMAP1 = t$umap1,
    L3UMAP2 = t$umap2
  )
  r[[cll]] <- t[[bestlds[c_j, "bestld"]]]
  rownames(r) <- r$barcode
  message(c_j, " has ", length(table(r[[cll]])), " cls.")
  return(r) 
}) |> do.call("rbind", args = _) |>
  x => `rownames<-`(x, x$barcode)

cellmeta <- merge(cellmetaL2, cellmetaL3, by = "barcode")
cellmeta$L1_2_3 <- with(cellmeta, L1_2 * 100 + L3)

# check # of L3
a <- map_int(rownames(bestlds), \(c_j) {
  fnm <- file.path(rdir, paste0(prefix, "_", c_j, ".cellmeta.leiden.csv"))
  t <- data.table::fread(fnm, sep = ",", header = TRUE, data.table = FALSE)
  bestld <- bestlds[c_j, "bestld"]
  cls <- t[[bestld]]
  r <- length(unique(cls))
  message(c_j, " has ", r, " cls.")
  return(r)
})
# 1041 clusters at L3 level

data.table::fwrite(
  x = cellmeta, file = file.path(rscdir,
    "pt.barcode2group.L3.20240401.csv"),
  row.names = FALSE, col.names = TRUE, sep = ","
)

# map L1_2_3 to noL4.
L123toL12 <- unique(cellmeta[ , c("L1_2_3", "L1_2")]) |>
  x => x[order(x$L1_2), ]
# L1_2 to noL4
L1_2noL4 <- map_int(seq_len(nrow(bestlds)), \(i) {
  c_j <- bestlds[i, "L1_2"]
  L1_2 <- as.integer(gsub("c", "", c_j))
  note <- bestlds[i, "note"]
  if (grepl("noL4", note)) {
    message(paste(c_j, L1_2, note))
    return(L1_2)
  }
  return(-5)
})
L1_2noL4 <- L1_2noL4[L1_2noL4 > 0]
L123noL4 <- map_int(seq_len(nrow(L123toL12)), \(i) {
  cj <- L123toL12[i, "L1_2_3"]
  ci <- L123toL12[i, "L1_2"]
  if (ci %in% L1_2noL4) {
    message(paste(ci, cj))
    return(cj)
  }
  return(-5)
})
L123noL4 <- L123noL4[L123noL4 > 0]
write.table(x = L123noL4,
  file = file.path(rscdir, "L3noL4.20240401.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE)

# * automatically sum L4
icll <- 4
pcll <- paste0("L", icll-1)
cll <- paste0("L", icll)
rdir <- file.path(workdir, "out", "tscc2", cll)
prefix <- "RNA_vg2000_npc30_k30_L1_2_3"
nmin <- 50

cellmetaL3 <- data.table::fread(
  file = file.path(rscdir, "pt.barcode2group.L3.20240401.csv"),
  sep = ",", header = TRUE, data.table = FALSE
) |> x => `rownames<-`(x, x$barcode)

L3noL4 <- read.table(file = file.path(rscdir, "L3noL4.20240401.txt"),
  header = F)$V1

c_js <- unique(cellmetaL3$L1_2_3)
n_cjs <- table(cellmetaL3$L1_2_3) |>
  as.data.frame(x = _, stringsAsFactors = FALSE) |>
  setNames(object = _, nm = c("L123", "n")) |>
  x => `rownames<-`(x, paste0("c", x$L123))
L3noL4_all <- union(L3noL4,
  as.integer(n_cjs$L123[n_cjs$n < 50]))
L4_cjs <- setdiff(c_js, L3noL4_all)

# get auto best leiden resos per L4
bestlds <- map_chr(L4_cjs, getBestlds) |>
  x => data.frame(
    L1_2_3 = paste0("c", L4_cjs),
    bestld = x,
    note="auto"
  ) |>
  x => x[order(x$L1_2_3), ]
  
data.table::fwrite(
  x = bestlds,
  file = file.path(rscdir, "auto.bestlds.L1_2_3.20240402.csv"),
  row.names = FALSE, col.names = TRUE
)

# summarize L4 after manual checking
raw_bestlds_sz <- data.table::fread(
  file = file.path(rscdir, "auto.bestlds.L1_2_3.20240405_SZ.csv"),
  header = TRUE, data.table = FALSE
)
needL5_sz <- map_lgl(seq_len(nrow(raw_bestlds_sz)), \(i) {
  if (grepl("noL5", raw_bestlds_sz[i, "note"])) {
    FALSE
  } else {
    TRUE
  }
})

raw_bestlds_zw <- data.table::fread(
  file.path(rscdir, "auto.bestlds.L1_2_3.20240402_ZW.csv"),
  header = TRUE, data.table = FALSE
)

# diff with sz: here optimal ones in the note
bestld_zw <- map_chr(seq_len(nrow(raw_bestlds_zw)), \(i) {
  if (grepl("auto", raw_bestlds_zw[i, "note"])) {
    return(raw_bestlds_zw[i, "bestld"])
  }
  return(str_split_1(raw_bestlds_zw[i, "note"], ";")[1])
})
needL5_zw <- map_lgl(seq_len(nrow(raw_bestlds_zw)), \(i) {
  if (grepl("noL5", raw_bestlds_zw[i, "note"])) {
    FALSE
  } else {
    TRUE
  }
})

bestld_sz <- data.frame(
  L1_2_3 = raw_bestlds_sz$L1_2_3,
  bestld = raw_bestlds_sz$bestld,
  needL5 = needL5_sz
)

bestld_zw <- data.frame(
  L1_2_3 = raw_bestlds_zw$L1_2_3,
  bestld = bestld_zw,
  needL5 = needL5_zw
)
bestlds <- rbind(bestld_sz, bestld_zw)
rownames(bestlds) <- bestlds$L1_2_3
data.table::fwrite(x =bestlds,
  file = file.path(rscdir, "auto.bestlds.L1_2_3.20240413.all.csv"),
  row.names = FALSE, col.names = TRUE, sep = ",")

# load leiden results on L4 then
cellmetaL4 <- lapply(rownames(bestlds), function(c_j) {
  fnm <- file.path(rdir, paste0(prefix, "_", c_j, ".cellmeta.leiden.csv"))
  t <- data.table::fread(fnm, sep = ",", header = TRUE, data.table = FALSE)
  bestld <- bestlds[c_j, "bestld"]
  if (! (bestld %in% colnames(t)) ) {
    stop(paste0(bestld, " not exist in ", c_j))
  }
  r <- data.frame(
    barcode = t$barcode,
    L4UMAP1 = t$umap1,
    L4UMAP2 = t$umap2
  )
  r[[cll]] <- t[[bestlds[c_j, "bestld"]]]
  rownames(r) <- r$barcode
  message(c_j, " has ", length(table(r[[cll]])), " cls.")
  return(r) 
}) |> do.call("rbind", args = _) |>
  x => `rownames<-`(x, x$barcode)

cellmetanoL4 <- data.frame(
  barcode = setdiff(rownames(cellmetaL3), rownames(cellmetaL4)),
  L4UMAP1 = NA,
  L4UMAP2 = NA,
  L4 = 0
)

cellmetaL4all <- rbind(cellmetanoL4, cellmetaL4)
rownames(cellmetaL4all) <- cellmetaL4all$barcode

cellmeta <- merge(cellmetaL3, cellmetaL4all, by = "barcode")
# length(unique(cellmeta[cellmeta$L1_2_3 == 90002, "L4"]))
cellmeta$L1_2_3_4 <- with(cellmeta, L1_2_3 * 100 + L4)
length(unique(cellmeta$L1_2_3_4))
data.table::fwrite(
  x = cellmeta, file = file.path(rscdir, "pt.barcode2group.L4.20240413.csv"),
  row.names = FALSE, col.names = TRUE, sep = ","
)

# * prepare L5-level clustering
# prepare L4 without needing for L5
# - only one cluster in L4
# - noL5 after manual check
# - *cluster with <= 50 cells (pipeline will consider this)
# - no L4 clustering by itself

cellmeta <- data.table::fread(
  file = file.path(rscdir, "pt.barcode2group.L4.20240413.csv"),
  sep = ",",
  header = TRUE,
  data.table = FALSE
) |> x => `rownames<-`(x, x$barcode)

# - no L4 clustering by itself
L4noL5_sinceL3 <- L3noL4_all * 100
# results from clustering on L1_2_3
L3toncl<- lapply(rownames(bestlds), function(c_j) {
  fnm <- file.path(rdir, paste0(prefix, "_", c_j, ".cellmeta.leiden.csv"))
  t <- data.table::fread(fnm, sep = ",", header = TRUE, data.table = FALSE)
  bestld <- bestlds[c_j, "bestld"]
  if (! (bestld %in% colnames(t)) ) {
    stop(paste0(bestld, " not exist in ", c_j))
  }
  ncl <- length(table(t[[bestlds[c_j, "bestld"]]]))
  data.frame(L123 = c_j, ncl = ncl)
}) |> do.call("rbind", args = _)

# - only one cluster in L4
# only 5 not in L4_noNeedL5
# due to manual annnot in note with "noL5"
L4noL5_onlyOneCluster <- with(L3toncl, L123[ncl < 2]) |>
  y => as.integer(gsub(pattern = "c", replacement = "", x = y)) * 100

L4toL3 <- unique(cellmeta[ , c("L1_2_3_4", "L1_2_3")]) |>
  x => x[order(x$L1_2_3), ] |>
  x => `rownames<-`(x, x$L1_2_3_4)
L3_noNeedL5 <- with(bestlds, L1_2_3[!needL5]) |>
  y => as.integer(gsub("c", "", x = y))

# - noL5 after manual check
# 2861
L4_noNeedL5 <- with(L4toL3, L1_2_3_4[L1_2_3 %in% L3_noNeedL5])

# 2935
L4noL5 <- union(x = as.integer(L4noL5_sinceL3),
  y = as.integer(L4noL5_onlyOneCluster)) |>
  union(x = _, y = L4_noNeedL5)

# - ncluster <= nmin
# - 145
# our pipeline will include this automatically
L4_nmin <- table(cellmeta$L1_2_3_4) |>
  y => dimnames(y)[[1]][y <= nmin]

write.table(x = L4noL5,
  file = file.path(rscdir, "L4noL5.20240413.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE)

# * summarize L5
icll <- 5
pcll <- paste0("L", icll-1)
cll <- paste0("L", icll)
rdir <- file.path(workdir, "out", "tscc2", cll)
prefix <- "RNA_vg2000_npc30_k30_L1_2_3_4"
nmin <- 50

cellmetaL4 <- data.table::fread(
  file = file.path(rscdir, "pt.barcode2group.L4.20240413.csv"),
  sep = ",", header = TRUE, data.table = FALSE
) |> x => `rownames<-`(x, x$barcode)

L4noL5 <- read.table(
  file = file.path(rscdir, "L4noL5.20240413.txt"),
  header = F
)$V1

c_js <- unique(cellmetaL4$L1_2_3_4)
n_cjs <- table(cellmetaL4$L1_2_3_4) |>
  as.data.frame(x = _, stringsAsFactors = FALSE) |>
  setNames(object = _, nm = c("L4", "n")) |>
  x => `rownames<-`(x, paste0("c", x$L4))
L4noL5_all <- union(L4noL5,
  as.integer(n_cjs$L4[n_cjs$n < nmin]))
L4_cjs <- setdiff(c_js, L4noL5_all)

# get auto best leiden resos per L4
bestlds <- map_chr(L4_cjs, getBestlds) |>
  x => data.frame(
    L1_2_3_4 = paste0("c", L4_cjs),
    bestld = x,
    note="auto"
  ) |>
  x => x[order(x$L1_2_3_4), ]
  
data.table::fwrite(
  x = bestlds,
  file = file.path(rscdir, "auto.bestlds.L1_2_3_4.20240415.csv"),
  row.names = FALSE, col.names = TRUE
)

# load manual check results
raw_bestlds <- data.table::fread(
  file = file.path(rscdir, "auto.bestlds.L1_2_3_4.20240415.csv"),
  sep = ",", header = TRUE, data.table = FALSE
)

tmp_bestld <- map_chr(seq_len(nrow(raw_bestlds)), \(i) {
  if (grepl("auto", raw_bestlds[i, "note"])) {
    return(raw_bestlds[i, "bestld"])
  }
  return(str_split_1(raw_bestlds[i, "note"], ";")[1])
})

bestlds <- data.frame(
  L4 = raw_bestlds$L1_2_3_4,
  bestld = tmp_bestld
)
rownames(bestlds) <- bestlds$L4

cellmetaL5 <- lapply(rownames(bestlds), function(c_j) {
  fnm <- file.path(rdir, paste0(prefix, "_", c_j, ".cellmeta.leiden.csv"))
  t <- data.table::fread(fnm, sep = ",", header = TRUE, data.table = FALSE)
  bestld <- bestlds[c_j, "bestld"]
  if (! (bestld %in% colnames(t)) ) {
    stop(paste0(bestld, " not exist in ", c_j))
  }
  r <- data.frame(
    barcode = t$barcode,
    L5UMAP1 = t$umap1,
    L5UMAP2 = t$umap2
  )
  r[[cll]] <- t[[bestlds[c_j, "bestld"]]]
  rownames(r) <- r$barcode
  message(c_j, " has ", length(table(r[[cll]])), " cls.")
  return(r) 
}) |> do.call("rbind", args = _) |>
  x => `rownames<-`(x, x$barcode)

cellmetanoL5 <- data.frame(
  barcode = setdiff(rownames(cellmetaL4), rownames(cellmetaL5)),
  L5UMAP1 = NA,
  L5UMAP2 = NA,
  L5 = 0
)
cellmetaL5all <- rbind(cellmetaL5,cellmetanoL5)
rownames(cellmetaL5all) <- cellmetaL5all$barcode
cellmeta <- merge(cellmetaL4, cellmetaL5all, by = "barcode")
# 4302 clusters
cellmeta$L1_2_3_4_5 <- with(cellmeta, L1_2_3_4 * 100 + L5)
data.table::fwrite(
  x = cellmeta, file = file.path(rscdir, "pt.barcode2group.L5.20240415.csv"),
  row.names = FALSE, col.names = TRUE, sep = ","
)

# * read meta data info for paper writing
# * 2025-07-10
ptcmfnm <- file.path(here::here(), "meta",
  "pairedtag.cell.meta.all.240626.csv")
rawptcm <- data.table::fread(ptcmfnm, header = T, sep = ",", data.table = F)
ptcm <- rawptcm[rawptcm$annotQuality == "Good", ]

nL1 <- table(rawptcm$cluster.l1.id)
nL1_ <- table(ptcm$cluster.l1.id)

# * check # of clusters for paper writing
# * 2026-01-14
ptcmfnm <- file.path(here::here(), "meta",
  "pairedtag.cell.meta.all.240626.csv")
rawptcm <- data.table::fread(ptcmfnm, header = T, sep = ",", data.table = F)
ptcm <- rawptcm[rawptcm$annotQuality == "Good", ]

ptcm[ptcm$isNeuL1, "pairedtagCluster"] |>
  unique() |>
  length()

ptcm[!ptcm$isNeuL1, "pairedtagCluster"] |>
  unique() |>
  length()

## table(rawptcm$annotQuality)

##                                 Good                    LQ_L5r_low_number       LQ_L5r_low_transferlabel_score 
##                              2588038                                18756                                  232 
## LQ_UMAP_disperse_at_class_per_region 
##                               115863 

