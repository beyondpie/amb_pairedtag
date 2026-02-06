suppressMessages({
  suppressWarnings({
    library(cicero)
    library(tidyverse)
    library(reticulate)
    library(future.apply)
    library(Matrix)
    library(GenomicRanges)
    library(rtracklayer)
    library(tmpRpkg)
  })
})

# * meta
plan(multicore, workers = 10)
projd <- here::here()
outd <- file.path(projd, "13.cicero", "out")
set.seed(2024)
mm10_genome <- read.table(
  file = file.path(projd, "meta", "mm10.chrom.sizes.lite"),
  header = F, sep = "\t"
)

# gft here is a GRanges
mm10_gtf <- rtracklayer::import(
  file.path(projd, "meta", "gencode.vM25.annotation.gtf")
)
# about 8% top pairs of conns selected then
cicero_model_cutoff <- 0.001

# * input and output
connsfnm <- file.path(outd, "conns.all.H3K27ac.csv")
connsout <- file.path(outd, "distal2proximal.conns.H3K27ac.0.001.bedpe")
d2gout <- file.path(outd, "CRE2gene.vM25.all.H3K27ac.0.001.csv")
RNAexpcsv <- file.path(
  projd, "data", "pairedtag_ann",
  "pt.RNAseq.pseudobulk.count.sc.allgenesymbol.pdDataFrame.csv"
)
RNAexprds <- file.path(
  projd, "data", "pairedtag_ann",
  "pt.RNAseq.pseudobulk.count.sc.allgenesymbol.rds"
)
H3K27acExph5 <- file.path(
  projd, "data", "pairedtag_ann",
  "pt.H3K27ac.pseudobulk.count.sc.bestSPM.h5"
)
H3K27acExprds <- file.path(
  projd, "data", "pairedtag_ann",
  "pt.H3K27ac.pseudobulk.count.sc.bestSPM.rds"
)
cordfnm <- file.path(
  outd, "all.H3K27ac.pdc.cor.csv")
shufcordfnm <- file.path(
  outd, "shuf.H3K27ac.pdc.cor.csv"
)

# * functions
transformPeakToGR <- function(peaks) {
  uPeaks <- unique(peaks)
  peak <- uPeaks |>
    strsplit(x = _, split = ":|-", perl = T) |>
    do.call(what = "rbind", args = _) |>
    tibble::as_tibble(x = _) |>
    setNames(object = _, c("chr", "startFrom", "endTo")) |>
    dplyr::mutate(
      startFrom = as.integer(startFrom),
      endTo = as.integer(endTo),
      name = uPeaks
    ) |>
    as.data.frame() |>
    x => `rownames<-`(x, x$name)
  return(peak)
}

# * main
conns <- data.table::fread(
  file = connsfnm, sep = ",", header = T, data.table = F
)

# filter by weight
conns <- subset(conns, coaccess >= cicero_model_cutoff)
peak1GR <- transformPeakToGR(conns$Peak1) |>
  tmpRpkg::loadGRfromBed(
    beds = _, colnms = c("chr", "start", "end", "name")
  )
peak2GR <- transformPeakToGR(conns$Peak2) |>
  tmpRpkg::loadGRfromBed(
    beds = _, colnms = c("chr", "start", "end", "name")
  )

# annot distal and proximal
promoters <- tmpRpkg::aroundGeneTSS(
  mm10_gtf,
  up = 1000, down = 1000,
  mcolnms = c("gene_id", "type", "gene_type", "gene_name")
)

peak1_proximal <- tmpRpkg::FindOvlpRegionInB(
  A = peak1GR, B = promoters
)
peak2_proximal <- tmpRpkg::FindOvlpRegionInB(
  A = peak2GR, B = promoters
)

conns$Peak1Annot <- "distal"
conns$Peak2Annot <- "distal"
conns$Peak1Annot[conns$Peak1 %in% peak1_proximal$name] <- "proximal"
conns$Peak2Annot[conns$Peak2 %in% peak2_proximal$name] <- "proximal"

# filter by distal-proximal pairs
conns_distal_and_proximal <- subset(
  conns,
  ((Peak1Annot == "distal") & (Peak2Annot == "proximal")) |
    ((Peak2Annot == "distal") & (Peak1Annot == "proximal"))
)

# unify as distal to proximal
conns_distal2proximal <- do.call(rbind,
  args =
    list(
      conns_distal_and_proximal[
        conns_distal_and_proximal$Peak1Annot == "distal",
      ],
      conns_distal_and_proximal[
        conns_distal_and_proximal$Peak2Annot == "distal",
      ] |>
        x => data.frame(
          Peak1 = x$Peak2,
          Peak2 = x$Peak1,
          coaccess = x$coaccess,
          Peak1Annot = x$Peak2Annot,
          Peak2Annot = x$Peak1Annot
        )
    )
) |>
  unique(x = _) |>
  # remove duplicated ones by selecting one coaccess score
  x => x[!duplicated(paste(x$Peak1, x$Peak2, sep = "_")), ]

# use promoter to replace proximal
# TODO: expland the proximal to genes (if multiple, list all of them)
## peak_proximal <- transformPeakToGR(conns_distal2proximal$Peak2) |>
##   loadGRfromBed(beds = _)
## gene_proximal <- findOverlaps2(
##   query = promoters, subject = peak_proximal)

# output
CREs_distal <- conns_distal2proximal$Peak1 |>
  strsplit(x = _, split = ":|-", perl = T) |>
  do.call(rbind, args = _)
CREs_proximal <- conns_distal2proximal$Peak2 |>
  strsplit(x = _, split = ":|-", perl = T) |>
  do.call(rbind, args = _)

CRE_names <- with(
  conns_distal2proximal,
  paste(Peak1, Peak2, sep = "->")
) |>
  x => data.frame(name = x)

bedpe <- cbind(
  CREs_distal, CREs_proximal, CRE_names,
  data.frame(score = conns_distal2proximal$coaccess)
)

data.table::fwrite(x = bedpe, file = connsout, sep = "\t", col.names = F)

# TODO: show cicero conns in subclass-sepecific manner

# * anlayze the ATAC and RNA co-expressions
# ** map proximal to genes
mapProximalToGene4CisConns <- function(conns_local) {
  group_by(.data = conns_local, Peak2) |>
    group_map(.data = _, ~ {
      proximals <- transformPeakToGR(peaks = .y$Peak2) |>
        loadGRfromBed(beds = _)
      genes <- tmpRpkg::FindOvlpRegionInB(promoters, proximals)
      gene <- tibble(
        gene_id = genes$gene_id,
        gene_name = genes$gene_name
      )
      distal <- tibble(
        CRE = .x$Peak1,
        score = .x$coaccess
      )
      r <- tidyr::expand_grid(distal, gene)
      return(r)
    }) |>
    do.call(what = rbind, args = _)
}

chroms <- strsplit(conns_distal2proximal$Peak1, ":|-", perl = T) |>
  vapply(X = _, \(i) {
    i[1]
  }, "chr10") |>
  unique()

d2gList <- future.apply::future_lapply(chroms, \(chr) {
  conns_local <- conns_distal2proximal[
    grepl(chr, conns_distal2proximal$Peak1),
  ]
  r <- mapProximalToGene4CisConns(conns_local)
  message(chr, " done.")
  return(r)
})
## NOTE: duplicated genes due to peaks map to the same genes
d2g <- do.call(rbind, d2gList) |>
  as.data.frame() |>
  x => x[, c("CRE", "gene_name", "gene_id", "score")]
data.table::fwrite(
  x = d2g,
  file = d2gout, sep = ",", col.names = T, row.names = F
)

# ** get pseudo-bulk RNA and DNA expression
# the data are generated in another python script
## RNA
tmpRpkg::setCondaEnv()
pd <- reticulate::import("pandas", convert = F)
pexp <- pd$read_csv(RNAexpcsv,
  sep = ",", header = 0L, index_col = 0L
)
pexpdf <- reticulate::py_to_r(pexp)
old_nms <- rownames(pexpdf)
nms <- gsub("annot.sc.", "", old_nms, perl = F) |>
  gsub(" |  ", "_", x = _, perl = F) |>
  gsub("-|/", "_", x = _, perl = F)
rownames(pexpdf) <- nms
saveRDS(pexpdf, RNAexprds)
## DNA
pH3K27acExp <- pd$read_hdf(H3K27acExph5,
  key = "H3K27ac", mode = "r"
)
pH3K27acExpdf <- reticulate::py_to_r(pH3K27acExp)
old_nms <- rownames(pH3K27acExpdf)
nms <- gsub("annot.sc.", "", old_nms, perl = F) |>
  gsub(" |  ", "_", x = _, perl = F) |>
  gsub("-|/", "_", x = _, perl = F)
rownames(pH3K27acExpdf) <- nms
saveRDS(pH3K27acExpdf, H3K27acExprds)

# ** prepare matrix for correlation
CRE2gene <- data.table::fread(
  file = d2gout, sep = ",", header = T, data.table = F
)

uCREs <- unique(CRE2gene$CRE)
uGenes <- unique(CRE2gene$gene_name)
genes <- uGenes[uGenes %in% colnames(pexpdf)]
CRE2gene <- CRE2gene[CRE2gene$gene_name %in% genes, ]

groups <- rownames(pH3K27acExpdf)
CREexpmat <- pH3K27acExpdf[groups, uCREs] |>
  as.matrix()
geneExpmat <- pexpdf[groups, genes] |>
  as.matrix()

# ** generate the shuffled background
## get gene.loci data.frame
geneLoci <- promoters |>
  x => x[match(genes, mcols(x)$gene_name)]
ndist.gene <- table(CRE2gene$gene_name)
## get shuffled CREs
shuffleCRE.gene <- furrr::future_map(names(ndist.gene), \(g) {
  d <- unique(CRE2gene[CRE2gene$gene_name == g, "CRE"])
  sample(setdiff(uCREs, d), ndist.gene[g], replace = F)
}, .options = furrr::furrr_options(seed = 1L), .progress = T)
shuffleCRE.gene <- unlist(shuffleCRE.gene)

## get random CREs
geneMatchShuffleCRE <- lapply(names(ndist.gene), \(g) {
  rep(g, ndist.gene[g])
}) |> unlist()

shufCRE2gene <- data.frame(
  CRE = shuffleCRE.gene,
  gene_name = geneMatchShuffleCRE
)

# ** fast pearson correlation
#' @param CREexpmat raw count group by CRE
#' @param geneExpmat raw count group by gene
#' geneExpmat: CPM matrix
getPearsonCor <- function(CRE2gene, CREexpmat, geneExpmat) {
  # make sure they share the same subclasses
  groupCRE <- rownames(CREexpmat)
  groupGene <- rownames(geneExpmat)
  if ((length(groupCRE) != length(groupGene)) ||
    (any(groupCRE != groupGene))) {
    message("rownames are not consistent for the two matrix.")
    message("align them.")
    groups <- intersect(groupCRE, groupGene)
    CREexpmat <- CREexpmat[groups, ]
    geneExpmat <- geneExpmat[groups, ]
  }
  message("Get CPM per group.")
  CRECPM <- (CREexpmat * 10^6) / rowSums(CREexpmat)
  geneCPM <- (geneExpmat * 10^6) / rowSums(geneExpmat)
  message("Perform nomalization")
  CREs <- CRE2gene[, 1]
  genes <- CRE2gene[, 2]
  cmat <- t(CRECPM[, CREs]) |>
    log1p(x = _) |>
    tmpRpkg::scaleByRow(mat = _)
  gmat <- t(geneCPM[, genes]) |>
    log1p(x = _) |>
    tmpRpkg::scaleByRow(mat = _)
  ## pearson
  message("Calculate Pearson Correlation.")
  tmpRpkg::fastCorOnRow.chunk(
    p1byc = cmat, p2byc = gmat, p12pairs = CRE2gene,
    corMethod = "pearson"
  )
}
cor_CRE2gene <- getPearsonCor(CRE2gene, CREexpmat, geneExpmat)
cor_shufCRE2gene <- getPearsonCor(shufCRE2gene, CREexpmat, geneExpmat)

cor2df <- function(r) {
  names <- names(r) |>
    strsplit(x = _, "@") |>
    do.call(rbind, args = _)
  data.frame(
    r1 = names[, 1],
    r2 = names[, 2],
    score = r
  )
}
cordf_CRE2gene <- cor2df(cor_CRE2gene) |> unique()
cordf_shufCRE2gene <- cor2df(cor_shufCRE2gene) |> unique()
data.table::fwrite(cordf_CRE2gene,
  file = cordfnm, col.names = F, row.names = F, sep = ",")
data.table::fwrite(cordf_shufCRE2gene,
  file = shufcordfnm, col.names = F, row.names = F, sep = ",")
