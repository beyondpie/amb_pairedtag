.libPaths("/opt/homebrew/Caskroom/miniforge/base/lib/R/library")
library(tidyverse)
library(ComplexHeatmap)

# * functions


loadsc2ddatCREsignal <- function(h) {
  d <- file.path(outd, str_glue("sc2distalDARChrACREovlpTE.{h}"))
  scs <- data.table::fread(
    file = file.path(d, "scs.bed"),
    header = F,
    data.table = F
  )$V1
  cres <- data.table::fread(
    file = file.path(d, "te.bed"),
    header = F,
    data.table = F,
    sep = "\t"
  ) |>
    setNames(object = _, nm = c("chrom", "startFrom", "endTo")) |>
    x => with(x,
      paste(chrom, paste(startFrom, endTo, sep="-"), sep = ":")) |>
  unlist()
  m <- data.table::fread(
    file = file.path(d, "mat.csv"),
    header = F, sep = ",", data.table = F
  ) |> as.matrix() |> t()
  colnames(m) <- cres
  rownames(m) <- scs
  m
}

genhm <- function(m, h, showsc = T,
                  colQmin = 0.05, colQmax = 0.99,
                  color = "red") {
  lowVal <- quantile(m, colQmin)
  highVal <- quantile(m, colQmax)
  # m <- mK27ac[ptscs, creord]
  col_fun <- circlize::colorRamp2(
    c(lowVal, highVal), c("white", color))
  legendLab <- c(round(lowVal, 1), round(highVal,1))

  Heatmap(m, name = "Distal DAR-ChrA CREs ovlp with TEs",
    col = col_fun,
    row_names_gp = gpar(fontsize = 5),
    cluster_rows = F,
    cluster_columns = F,
    use_raster = T,
    show_row_names = showsc,
    row_names_side = "left",
    show_column_names = F,
    heatmap_legend_param = list(
      title = str_glue("logRPKM {h}"),
      at = c(lowVal, highVal),
      labels = legendLab,
      direction = "horizontal"
    )
  )
}

# * main
projd <- here::here()
workd <- file.path(projd, "17.repressiveMarks")
outd <- file.path(workd, "out", "ovlpTE")

allCREfnm <- file.path(projd, "data", "CRE",
  "allCRE.amb.PairedTag.annot.tsv")

allCREs <- data.table::fread(
  file = allCREfnm, header = T, sep = "\t", data.table = F)

# dda: distal, DAR, Chr-A
index <- with(allCREs, isDAR == "DAR" &
                         isDistal == "distal" &
                         chromHMMState == "Chr-A")
# 470,164
ddaCREs <- allCREs[index, ] |>
  x => with(x,
    paste(chrom, paste(startFrom, endTo, sep = "-"), sep = ":")) |>
  unique()

CREovlpTEfnm <- file.path(workd, "out", "ovlpTE",
  "CRE-TE-TECat.tsv")
CREovlpTE <- data.table::fread(
  file = CREovlpTEfnm, sep = "\t", header = F, data.table = F
) |>
  x => paste(x[,1], paste(x[,2], x[,3], sep = "-"), sep = ":") |>
  unique()

# ddat: distal, DAR, Chr-A, ovlp with TE
# 105,717
ddatCREs <- intersect(ddaCREs, CREovlpTE) |>
  sort()
data.table::fwrite(x = data.frame(cre = ddatCREs),
  file = file.path(outd, "ovlpTE", "distalDARChrACREs.ovlpwithTE.txt"),
  row.names = F, col.names = F
  )

# * save subclass 2 ddat CRE list
# for each subclass
scs <- allCREs[index, "subclass"] |>
  unique() |>
  sort()
invisible(
  lapply(scs, \(sc) {
    fnm <- file.path(workd, "out", "sc2ddatCRE",
      str_glue("{sc}_distalDARChrA-ovlpTE_CRE.bed"))
    t <- allCREs[index, ] |>
      x => x[x$subclass == sc, 1:3] |>
      x => setNames(object = x, nm = c("chrom", "startFrom", "endTo"))
    t$chromNum <- t[ , 1] |>
      x => gsub("chr", "", x) |>
      x => as.integer(x)
    t <- t[order(t$chromNum, t$startFrom, t$endTo),
      c("chrom", "startFrom", "endTo")]
    
    tstr <- with(t, paste(chrom,
      paste(startFrom, endTo, sep = "-"), sep = ":"))
    t <- t[which(tstr %in% ddatCREs), ]
    if (nrow(t) >= 1) { 
      data.table::fwrite(
        x = t, file = fnm, row.names = F, col.names = F, sep = "\t")
    } else {
      message(str_glue("{sc} has no CREs left."))
    }
  })
)

# * get subclass 2 ddatCREs matrix with
# different histone modifications
mK27ac <- loadsc2ddatCREsignal(h = "H3K27ac")

# set the column orders for ddatCREs
ptscs <- rownames(mK27ac)
sc_creord <- lapply(ptscs, \(sc) {
  indexsc <- allCREs$subclass == sc
  x <- allCREs[index & indexsc, ] |>
    x => with(x, paste(chrom,
      paste(startFrom, endTo, sep = "-"), sep = ":")) |>
    x => intersect(x, colnames(mK27ac))
  s <- mK27ac[sc, x] |>
    y => order(y, decreasing = T)
  data.frame(
    sc = sc,
    cre = x[s]
  )
}) |>
  do.call(what = rbind, args = _)

ntop <- 20
creord <- c()
for (sc in ptscs) {
  cur <- creord
  x <- setdiff(sc_creord[sc_creord$sc == sc, "cre"], creord)
  if (length(x) == 1) {
    if (!(x %in% cur)) {
      creord <- c(cur, x)
    }
  }
  if (length(x) > 1) {
    y <- x[1: min(length(x), ntop)] |>
      setdiff(x = _, y = creord)
    if (length(y) > 0) {
      creord <- c(cur, y)
    }
  }
}


# ** load sc2CRE histone activities
mK27me3 <- loadsc2ddatCREsignal(h = "H3K27me3")
mK9me3 <- loadsc2ddatCREsignal(h = "H3K9me3")
mK4me1 <- loadsc2ddatCREsignal(h = "H3K4me1")

saveRDS(
  object = list(
    ptscs = ptscs,
    cre4plot = creord,
    mK27ac = mK27ac,
    mK27me3 = mK27me3,
    mK4me1 = mK4me1,
    mK9me3 = mK9me3
  ), file = file.path(outd, "saveHeatmapData4sc2DDARCREovlpTE.rds"))

# ** plot
hmK27ac <- genhm(
  m = log1p(mK27ac[ptscs, creord]),
  h = "H3K27ac",
  showsc = T,
  colQmin = 0.05, colQmax = 0.99,
  color = "#D85641"
  )

hmK27me3 <- genhm(
  m = log1p(mK27me3[ptscs, creord]),
  h = "H3K27me3",
  showsc = T,
  colQmin = 0.05, colQmax = 0.99,
  color = "#E7AA23"
  )

hmK4me1 <- genhm(
  m = log1p(mK4me1[ptscs, creord]),
  h = "H3K4me1",
  showsc = T,
  colQmin = 0.05, colQmax = 0.99,
  color = "#90C03E"
  )

hmK9me3 <- genhm(
  m = log1p(mK9me3[ptscs, creord]),
  h = "H3K9me3",
  showsc = T,
  colQmin = 0.05, colQmax = 0.99,
  color = "#5551A2"
  )

hm <- hmK27ac + hmK4me1 + hmK27me3 + hmK9me3

withr::with_pdf(
  new = file.path(outd, str_glue("heatmap_sc2ddatCRE_all.pdf")),
  code = {draw(hm)},
  width = 30,
  height = 10
)
