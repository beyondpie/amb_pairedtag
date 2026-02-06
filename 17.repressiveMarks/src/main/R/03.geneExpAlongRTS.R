library(tidyverse)
library(tmpRpkg)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

# * meta
projd <- here::here()
workd <- file.path(projd, "17.repressiveMarks")
rtsd <- file.path(workd, "out")
outd <- file.path(workd, "out")
figd <- file.path(workd, "figure")

# mod <- "H3K27me3"
mod <- "H3K9me3"
qtop <- 0.05

# load Allen's Gene lists
AllenGened <- file.path(projd, "meta", "wmbAllenUpdate")
tfDE <- data.table::fread(file = file.path(AllenGened, "TF534.txt"),
  head = F, data.table = F)$V1
tfMod <- data.table::fread(file = file.path(AllenGened,
  "TF261CoExpModule52.tsv"),
  head = F, sep = "\t", data.table = F) |>
  setNames(object = _, nm = c("mod", "gene", "score")) |>
  x => `rownames<-`(x, x$gene)

# * functions
loadRTS <- function(fnm) {
  data.table::fread(
    file = fnm, sep = "\t", header = F, data.table = F
  ) |> setNames(object = _,
    nm = c("chrom", "startFrom", "endTo", "gene", "rts", "strand"))
}

getAvgExpAlongBin <- function(xscbyg, genes, binSize = 100) {
  breaks <- seq(from = 1, to = length(genes), by = binSize)
  vapply(seq_len(length(breaks)-1), \(i) {
    g <- intersect(genes[breaks[i] : breaks[i+1]], colnames(xscbyg))
    if(length(g) < 1) {
      NULL
    } else {
      mean(xscbyg[ , g])
    }
  }, 0.0)
}

getGeneEnrichRatioAlongBin <- function(gRTS, gTF, binSize = 200) {
  breaks <- seq(from = 1, to = length(gRTS), by = binSize)
  vapply(seq_len(length(breaks)-1), \(i) {
    g <- intersect(gRTS[breaks[i]:breaks[i+1]], gTF)
    length(g) / length(gTF)
  }, 0.0)
}

# * main
# load RTS
rts <- loadRTS(file.path(rtsd,
  str_glue("gene2RTS.{mod}.top{qtop}.all.bed"))) |>
  x => x[!duplicated(x$gene), ] |>
  x => `rownames<-`(x, x$gene)

# load pseudo-bulk level gene expressions
logcpmscbyg <- readRDS(file.path(projd, "data", "pairedtag_ann",
  "pt.RNAseq.pseudobulk.CPM.scbyg.rds")) |>
  x => log1p(x)

avgExpAlongBin <- getAvgExpAlongBin(logcpmscbyg, genes = rts$gene, binSize = 100)

plot(x = seq_along(avgExpAlongBin)[1:120], y = avgExpAlongBin[1:120])

# 0.78 for H3K27me3, -0.34 for H3K9me3
cor(x = seq_along(avgExpAlongBin)[1:120],
  y = avgExpAlongBin[1:120], method = "spearman") 

# TF enrichment along the bins
tfDE <- intersect(tfDE, rts$gene)
tfDEEnrichRatio <- getGeneEnrichRatioAlongBin(rts$gene, tfDE, binSize  = 200)
plot(x = seq_along(tfDEEnrichRatio), y = tfDEEnrichRatio)

tfModGene <- intersect(tfMod$gene, rts$gene)
tfModRatio <- getGeneEnrichRatioAlongBin(rts$gene, tfModGene, binSize  = 200)
plot(x = seq_along(tfModRatio), y = tfModRatio)




