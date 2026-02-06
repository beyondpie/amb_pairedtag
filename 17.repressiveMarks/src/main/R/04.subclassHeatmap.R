library(tidyverse)
library(tmpRpkg)
library(ComplexHeatmap)
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

# * functions
loadRTS <- function(fnm) {
  data.table::fread(
    file = fnm, sep = "\t", header = F, data.table = F
  ) |> setNames(object = _,
    nm = c("chrom", "startFrom", "endTo", "gene", "rts", "strand"))
}

# * main
# load RTS
rts <- loadRTS(file.path(rtsd,
  str_glue("gene2RTS.{mod}.top{qtop}.all.bed"))) |>
  x => x[!duplicated(x$gene), ] |>
  x => `rownames<-`(x, x$gene)

# get topRTS
dy <- abs(diff(rts$rts))
cumdy <- cumsum(dy)
ng <- min(which(cumdy >= 0.9))
minRTS <- rts$rts[ng]
topRTS <- rts[rts$rts >= minRTS, ]

# load pseudo-bulk level gene expressions
logcpmscbyg <- readRDS(file.path(projd, "data", "pairedtag_ann",
  "pt.RNAseq.pseudobulk.CPM.scbyg.rds")) |>
  x => log1p(x)

# prepare matrix for heatmap
geneTopRTS <- intersect(topRTS$gene, colnames(logcpmscbyg))
X <- logcpmscbyg[ , geneTopRTS] |>
  scale(x = _, center = T, scale = T)

scs <- rownames(logcpmscbyg)
combscs <- combn(scs, 2)

cors <- vapply(seq_len(ncol(combscs)), \(i) {
  cor(X[combscs[1, i], ], X[combscs[2, i], ], method = "pearson")
}, 0.0)

hmX <- matrix(data = 0.0, nrow = length(scs), ncol = length(scs),
  dimnames = list(
    scs, scs
  )
)
for (i in seq_len(ncol(combscs))) {
  s1 <- combscs[1, i]
  s2 <- combscs[2, i]
  hmX[s1, s2] <- cors[i]
  hmX[s2, s1] <- cors[i]
}

low.val.col <- 0.0
high.val.col <- 0.7
hmCor <- tmpRpkg::setHeatmap(
  mat = hmX,
  low.val.col = low.val.col,
  high.val.col = high.val.col ,
  col_fun = circlize::colorRamp2(seq(low.val.col, high.val.col, length = 60),
    viridis::viridis(60)
    ),
  showRowNames = F,
  showColNames = F,
  cluster_columns = F,
  cluster_rows = F,
  use_raster = T,
  legend_title = "PearsonCor with topRTS"
)

hmCor

