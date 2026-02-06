# 3. motif enrichment directly from gimme score, which is logoddratio
# - load motif-TF-subclass-motifscore
# - for each subclass, generate the TF motif-enrichment based score
# - heatmap

suppressMessages({
  suppressWarnings({
    library(tidyverse)
    # detach("package:tmpRpkg", unload = T)
    library(tmpRpkg)
    library(ComplexHeatmap)
    library(Seurat)
    library(future)
    library(future.apply)
  })
})

options(future.globals.maxSize = 50 * 1024^3)
future::plan(multicore, workers = 30)
# * meta
projd <- tmpRpkg:::projd
workd <- file.path(projd, "16.celloracle")
eRegulond <- file.path(workd, "out", "eRegulon")
outd <- file.path(workd, "out")
ptscMetafnm <- file.path(
  projd, "meta", "PairedTagSubclassMetaFromCellMetaChromHMM.csv"
)
ptscMeta <- data.table::fread(
  file = ptscMetafnm,
  sep = ",", header = TRUE, data.table = FALSE
) |>
  x => `rownames<-`(x, x$PairedTagName) |>
  x => x[x$ATAC > 0, ]

ptscs <- ptscMeta$PairedTagName
atacscs <- ptscMeta$ATACName |>
  x => gsub("^\\d+_", "", x)
sc_pt2atac <- data.frame(
  pt = ptscs,
  atac = atacscs
) |> x => `rownames<-`(x, x$pt)

# get number of CRE per subclass
ATACREd <- file.path(projd, "data", "snATAC", "subclass2CRE")
sc2nCRE <- lapply(ptscs, \(ptsc) {
  atacsc <- sc_pt2atac[ptsc, "atac"]
  fnm <- file.path(ATACREd, str_glue("{atacsc}.cCREs.txt"))
  nrow(
    data.table::fread(file = fnm, header = F, sep = ",", data.table = F)
  )
}) |>
  unlist() |>
  setNames(object = _, nm = ptscs)

# get subclass2tf from gene expression
figd <- file.path(projd, "99.figures/out/fig4")
sc2tfRNAmat <- readRDS(file.path(figd,
  "panelE.subclassSpecEnhTF.pairedTagRNA.scaledlogCPM.rds"))

# * functions
loadSubclassRegulon <- function(sc) {
  fnm <- file.path(eRegulond, str_glue("gimmeMotif.{sc}.CRE2tf.csv"))
  data.table::fread(
    file = fnm,
    sep = ",",
    header = T,
    data.table = F
  )
}

getSubclassEnrichScore <- function(eRegulon, sc) {
  eRegulon |>
    group_by(tf) |>
    summarise(
      subclass = sc,
      sumGimmeScore = sum(score),
      maxGimmeScore = max(score),
      avgGimmeScore = mean(score),
      motiftimes = sum(score > 0),
      sumATACPM = sum(atacCPM),
      mostEnrichMotif = motif_id[which.max(score)]
    ) |>
    as.data.frame()
}

transformMotifEnrich2mat <- function(motifEnrich, whichScore) {
  tmp <- with(motifEnrich, {
    data.frame(
      x = subclass,
      y = tf,
      z = motifEnrich[[whichScore]]
    )
  })
  with(tmp, {
    xnms <- unique(x)
    ynms <- unique(y)
    out <- matrix(
      data = 0.0, nrow = length(xnms), ncol = length(ynms),
      dimnames = list(xnms, ynms)
    )
    out[cbind(x, y)] <- z
    out
  })
}

normalizeMotifEnrich <- function(mat_scBytf, sc2nCRE) {
  scs <- rownames(mat_scBytf)
  diag(1.0 / sc2nCRE[scs]) %*% mat_scBytf |>
    log1p() |>
    x => `rownames<-`(x, scs)
}

# * main
# 1. get all the motif enrichment record
sc2mtfenrm <- future_lapply(ptscs, \(sc) {
  e <- loadSubclassRegulon(sc)
  s <- getSubclassEnrichScore(e, sc)
  message(sc, " done.")
  return(s)
}) |>
  setNames(object = _, nm = ptscs)
saveRDS(sc2mtfenrm, file.path(outd, "all.motifenrich.sum.rds"))

# 2. merge into one unifed data.frame
motifEnrich <- do.call(rbind, sc2mtfenrm)

# 3. get subclass2tf enrichment score matrix
m <- transformMotifEnrich2mat(
  motifEnrich,
  "maxGimmeScore"
) |>
  x => x[rownames(sc2tfRNAmat), colnames(sc2tfRNAmat)]


m <- transformMotifEnrich2mat(
  motifEnrich,
  "sumGimmeScore"
) |>
  x => x[rownames(sc2tfRNAmat), colnames(sc2tfRNAmat)] |>
  x => normalizeMotifEnrich(x, sc2nCRE)

up <- quantile(m, 0.90)
m[ m > up] <- up

# 4. drawHeatmap
low.val.col <- quantile(m, 0.1, na.rm = T)
high.val.col <- quantile(m, 0.9, na.rm = T)

legend_labels <- c(round(low.val.col, 1),
  round(high.val.col, 1))

col_fun <- circlize::colorRamp2(
  seq(low.val.col, high.val.col, length = 30),
  viridis::viridis(30)
)

p <- ComplexHeatmap::Heatmap(
  matrix = m,
  col = col_fun,
  cluster_columns = F,
  cluster_rows = F,
  show_row_names = T,
  row_names_gp = grid::gpar(fontsize = 7),
  column_names_gp = grid::gpar(fontsize = 7),
  show_column_names = T,
  top_annotation = NULL,
  left_annotation = NULL,
  use_raster = T,
  show_heatmap_legend = T,
  heatmap_legend_param = list(
    title = "motif Enrichment Scoare",
    at = c(low.val.col, high.val.col),
    labels = legend_labels,
    direction = "horizontal"
  )
)

withr::with_pdf(
  new = file.path(figd,
    "panelE.subclassSpecTFMotifEnrich.pdf"),
  code = {
    print(p)
  },
  width = 10,
  height = 12
)
