library(tidyverse)
library(future.apply)
future::plan(multisession, workers = 4)

# meta
projd <- "/mnt/tscc2/szu/projects/amb_pairedtag/"
setwd(projd)
workd <- file.path(projd, "04.peakcalling")
rscd <- file.path(workd, "src/main/resource")
peakd <- file.path(workd, "out", "macs3_nolambda")
mods <- c("H3K27ac", "H3K27me3", "H3K4me1", "H3K9me3")

# load metadata
neuscStat <- data.table::fread(
  file = file.path(rscd, "Neu.subclass.SexRep.stat.csv"),
  header = TRUE,
  sep = ",",
  data.table = FALSE
)
colnames(neuscStat) <- c(
  "group", "modality", "sexrep", "ncell", "nreads", "ratio")

nnspStat <- data.table::fread(
  file = file.path(rscd, "NN.supertype.SexRep.stat.csv"),
  sep = ",",
  header = TRUE,
  data.table = FALSE
)
colnames(nnspStat) <- c(
  "group", "modality", "sexrep", "ncell", "nreads"
)

groupStat <- rbind(neuscStat[ , 1:5], nnspStat)
t <- vapply(groupStat$group, \(v) {
  gsub(" +", "_",v) |>
    x => gsub("/", "_", x) |>
    x => gsub("-", "_", x)
}, "i")
rownames(groupStat) <- paste(t, groupStat$modality, 
                             groupStat$sexrep, sep = "-")


# functions
get_num_peaks <- function(group, mod, sex, qvalue = 0.01) {
  message(str_glue("{group}-{mod}-{sex}"))
  suffix <- ifelse(mod == "H3K27ac", "narrowPeak", "broadPeak")
  prefix <- str_glue("{group}-{mod}-{sex}")
  names <- c(
    #str_glue("{group}-{mod}-{sex}_peaks.{suffix}"),
    str_glue("{group}-{mod}-{sex}A_peaks.{suffix}"),
    str_glue("{group}-{mod}-{sex}B_peaks.{suffix}")
   # str_glue("{group}-{mod}-{sex}-shufA_peaks.{suffix}"),
   # str_glue("{group}-{mod}-{sex}-shufB_peaks.{suffix}")
   )
  fnms <- vapply(names, \(f) {file.path(peakd, group, f)}, "")
  npeaks <- vapply(fnms, \(f) {
    if (!file.exists(f)) {
      return(0)
    }
    if (file.size(f) < 5) {
      return(0)
    }
    d <- data.table::fread(file = f, header = FALSE, sep = "\t", 
                           data.table = FALSE)
    d <- d[d[,9] >= -log10(qvalue), ]
    nrow(d)
  }, 1)
  ncell <- vapply(c("A", "B"), \(i) {
    nm <- paste0(prefix, i)
    if (nm %in% rownames(groupStat)) {
      return(groupStat[nm, "ncell"])
    }
    0}, 0)
  nreads <- vapply(c("A", "B"), \(i){
    nm <- paste0(prefix, i)
    if (nm %in% rownames(groupStat)) {
      return(groupStat[nm, "nreads"])
    }
    0}, 0)
  data.frame(
    names = paste0(prefix, c("A", "B")),
    npeak = npeaks,
    ncell = ncell,
    nreads = nreads
  )
}

# main
groups <- list.files(peakd, no.. = TRUE)
# H3K27me3_Male <- future_lapply(groups, \(g) {
#   get_num_peaks(group = g, mod = "H3K27me3", sex = "Male")}) |> 
#   x => do.call(rbind, x)

peakCallingStat <- lapply(mods, \(m) {
  lapply(c("Male", "Female"), \(s) {
    future_lapply(groups, \(g) {
      get_num_peaks(group = g, mod = m, sex = s)
    }) |> x => do.call(rbind, x)
  }) |> x => do.call(rbind, x)
}) |> x => do.call(rbind, x)

data.table::fwrite(
  x = peakCallingStat,
  file = file.path(rscd, "macs3PeakCallingStat.csv"),
  sep = ",",
  col.names = TRUE,
  row.names = FALSE
)

plotRelation <- function(mod, x = "nreads", y = "npeak") {
  t <- peakCallingStat[grepl(mod, peakCallingStat$names), c(x, y)]
  colnames(t) <- c("X", "Y")
  t <- t[t$X > 0, ]
  t$X <- log10(t$X+1)
  t$Y <- log10(t$Y+1)
  ggplot(data = t, aes(x = X, y = Y)) +
    geom_point() +
    xlab(str_glue("log10 {x}")) +
    ylab(str_glue("log10 {y}")) +
    ggtitle(mod)
}
(p_H3K27ac <- plotRelation("H3K27ac", "nreads", "npeak"))
(p_H3K27me3 <- plotRelation("H3K27me3", "nreads", "npeak"))
(p_H3K4me1 <- plotRelation("H3K4me1", "nreads", "npeak"))
(p_H3K9me3 <- plotRelation("H3K9me3", "nreads", "npeak"))


# for debug
# get_num_peaks(group = "116_COAa_PAA_MEA_Barhl2_Glut",
#              mod = "H3K27me3", sex = "Male")

# add corresponding cell number
# filter null data.frame (check why they are null)
