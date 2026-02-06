library(tidyverse)

projd <- here::here()
workd <- file.path(projd, "17.repressiveMarks")
fromKanglid <- file.path(workd, "out", "fromKangli")

# load pdc / ppdc
# all CREs
peakAnno <- readRDS(file.path(fromKanglid, "peak.anno.rename.rds"))

# peak Ovlp with TEs
peakAnnoTE <- readRDS(file.path(fromKanglid, "peak.anno.TEs.rds"))

# link TE with CRE (Fig5.f)
# some genes are NA, which should be from the results of
# joint two different tables
# This is for highTESubclass-enriched TEs
ppdcTE2CRE <- readRDS(file.path(fromKanglid,
  "TE.cCREs_pos.pdc_highTE_vs_other_test.rds"))

# TE family, superfamily, and detailed annotation
TEdict <- readRDS(file.path(fromKanglid, "TE_dict.rds"))


# load subclass-level ChrA-CREs
ptCREs <- data.table::fread(
  file.path(projd, "data", "CRE", "allCRE.amb.PairedTag.annot.tsv"),
  header = T, sep = "\t", data.table = F)

highTEscs <- data.table::fread(
  file.path(projd, "meta", "highTESubclass.csv"),
  header = T, sep = ",", data.table = F
) |>
  x => x[, 2] |>
  x => unlist(x)

highTEChrACREs <- ptCREs[
  (ptCREs$subclass %in% highTEscs) & (ptCREs$chromHMMState == "Chr-A"), ] |>
  x => with(x, paste(chrom, paste(startFrom, endTo, sep = "-"), sep = ":")) |>
  x => unique(x)
  
# get highTE-ChrACRE ovlp TEs
CRE2TE <- data.table::fread(
  file.path(workd, "out", "ovlpTE", "CRE-TE-TECat.tsv"),
  header = F, sep = "\t", data.table = F
) |>
  setNames(object = _,
    nm = c("CRE_chrom", "CRE_startFrom", "CRE_endTo",
      "TE_chrom", "TE_startFrom", "TE_endTo",
      "TEFamily", "TEHomerFam", "TESubfamily"
      )) |>
  mutate(CREStr = paste(CRE_chrom,
    paste(CRE_startFrom, CRE_endTo, sep = "-"), sep = ":"))

# This is CRE overlapped with TE
highTEChrATE <- CRE2TE[
  CRE2TE$CREStr %in% highTEChrACREs, "CREStr"] |>
  x => unique(x)
  
# filter ppdcTE2CRE
# activate TE to genes
at2g <- ppdcTE2CRE[ppdcTE2CRE$TEs %in% highTEChrATE, ]
at2g$DCA <- "Non-DCA"
at2g$DCA[
  (at2g$TE.cCREs_pos.pdc_highTE_vs_other_allTestP_fdr <= 0.01) &
  (at2g$TE.cCREs_pos.pdc_highTE_vs_other_log2OR >= 1.5)] <- "DCA"

at2g$neglog10FDR <- (-log10(at2g$TE.cCREs_pos.pdc_highTE_vs_other_allTestP_fdr))
at2g$log2fd <- at2g$TE.cCREs_pos.pdc_highTE_vs_other_log2OR

# filtered activate TE to genes by logfd and fpr
## fat2g <- at2g[
## (at2g$TE.cCREs_pos.pdc_highTE_vs_other_allTestP_fdr <= 0.01) &
##   (at2g$TE.cCREs_pos.pdc_highTE_vs_other_log2OR >= 1.5) &
##   (!is.na(at2g$gene)),
## ]

fat2g <- at2g[
(at2g$TE.cCREs_pos.pdc_highTE_vs_other_allTestP_fdr <= 0.01) &
  (at2g$TE.cCREs_pos.pdc_highTE_vs_other_log2OR >= 1.5), 
]

saveRDS(
  object = at2g,
  file = file.path(workd, "out/drawFig5", "HighTE_ChrATE.vsOtherCREs.rds"))

examplePools <- fat2g[!is.na(fat2g$gene), ]
saveRDS(
  object = examplePools,
  file = file.path(workd, "out/drawFig5", "TE2Gene.forExamples.rds")
)

# check this TE2Gene examples.
r <- readRDS(file.path(workd, "out/drawFig5", "TE2Gene.forExamples.rds"))


# * draw Fig5f
DCAcolors <- c("gray", "red") |>
  setNames(object = _, nm = c("Non-DCA", "DCA"))
TEShapes <- c(24, 21, 23, 22) |>
  setNames(object = _, nm = c("DNA", "LINE", "LTR", "SINE"))

(
  p <- ggplot(data = at2g, aes(x = log2fd, y = neglog10FDR,
    fill = DCA, color = DCA, shape = AnnoSuperFam)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_hline(yintercept = 2, color = "black", linetype="dashed") +
    geom_vline(xintercept = 1.5, color = "black", linetype = "dashed")+
    scale_fill_manual(values = DCAcolors) +
    scale_color_manual(values = DCAcolors) +
    scale_shape_manual(values = TEShapes) +
    theme_bw(base_size = 25)
)

# * Fig5g
# for HOMER enrichment analysis
t_chr2coords <- lapply(fat2g$TEs, \(g) {
  r = str_split_1(g, pattern = ":|-")
  data.frame(
    chr = r[1],
    s = as.integer(r[2]),
    e = as.integer(r[3])
  )
}) |>
  do.call(what = rbind, args = _) |>
  mutate(
    name = fat2g$AnnoDetail,
    score = 0.0,
    strand = "+"
  ) |>
  unique()

data.table::fwrite(
  x = t_chr2coords,
  file = file.path(workd, "out/drawFig5", "teCRE4homer.bed"),
  col.names = F,
  row.names = F,
  sep = "\t"
)

bg <- lapply(at2g$TEs[at2g$DCA == "Non-DCA"], \(g) {
  r = str_split_1(g, pattern = ":|-")
  data.frame(
    chr = r[1],
    s = as.integer(r[2]),
    e = as.integer(r[3])
  )
}) |>
  do.call(what = rbind, args = _) |>
  mutate(
    name = at2g$AnnoDetail[at2g$DCA == "Non-DCA"],
    score = 0.0,
    strand = "+"
  ) |>
  unique()

data.table::fwrite(
  x = bg,
  file = file.path(workd, "out/drawFig5", "bg_teCRE4homer.bed"),
  col.names = F,
  row.names = F,
  sep = "\t"
)
  

# * get subclass-specific ChrA-TEs
# 2026-01-20
projd <- here::here()
workd <- file.path(projd, "17.repressiveMarks")
fromKanglid <- file.path(workd, "out", "fromKangli")

ptCREs <- data.table::fread(
  file.path(projd, "data", "CRE",
    "allCRE.amb.PairedTag.annot.tsv"),
  header = T, sep = "\t", data.table = F)

DARChrACREs <- ptCREs[
  (ptCREs$chromHMMState == "Chr-A") &
  (ptCREs$isDistal == "distal") &
  (ptCREs$isDAR == "DAR"), ] |>
  x => mutate(.data = x, coord = paste(x$chrom,
    paste(x$startFrom, x$endTo, sep = "-"), sep = ":")) 
  

CRE2TE <- data.table::fread(
  file.path(workd, "out", "ovlpTE", "CRE-TE-TECat.tsv"),
  header = F, sep = "\t", data.table = F
) |>
  setNames(object = _,
    nm = c("CRE_chrom", "CRE_startFrom", "CRE_endTo",
      "TE_chrom", "TE_startFrom", "TE_endTo",
      "TEFamily", "TEHomerFam", "TESubfamily"
      )) |>
  mutate(CREStr = paste(CRE_chrom,
    paste(CRE_startFrom, CRE_endTo, sep = "-"), sep = ":"))

ptsc <- unique(DARChrACREs$subclass)
scDARChrATEs <- lapply(ptsc, \(sc) {
  s <- DARChrACREs[DARChrACREs$subclass == sc, "coord"]
  r <- CRE2TE[CRE2TE$CREStr %in% s, ]
  data.frame(
    TE = with(r, paste(TE_chrom,
      paste(TE_startFrom, TE_endTo, sep = "-"), sep = ":")),
    cCRE = r$CREStr,
    subclass = sc,
    TEFamily = r$TEFamily,
    TESubfamily = r$TESubfamily
  )
}) |>
  do.call(what = rbind, args = _)
data.table::fwrite(x = scDARChrATEs,
  file = file.path(workd, "out", "ambpt_subclassDARChrATEs.csv"),
  col.names = T, row.names = F, sep = ","
  )

# * get subclass-specific ChrA-CREs
scDARChrACREs <- DARChrACREs[ , c("coord", "subclass")] |>
  unique()
colnames(scDARChrACREs) <- c("cCRE", "subclass")

data.table::fwrite(x = scDARChrACREs,
  file = file.path(workd, "out", "ambpt_subclassDARChrACREs.csv"),
  col.names = T, row.names = F, sep = ","
  )

dsCREbws <- readRDS(file.path(projd, "06.ChromHMM",
  "out", "PanelD.ds1000.rds"))

scEnhs <- lapply(seq_along(colnames(dsCREbws$K27ac)), \(i) {
  atacsc <- colnames(dsCREbws$ATAC)[i]
  ptsc <- colnames(dsCREbws$K27ac)[i]
  s <- scDARChrACREs[scDARChrACREs$subclass == ptsc, ]
  cres <- intersect(s$cCRE, rownames(dsCREbws$ATAC))
  atac <- dsCREbws$ATAC[cres, atacsc]
  k27 <- dsCREbws$K27ac[cres, ptsc]
  ind <- order(atac, decreasing = T)
  data.frame(
    cCRE = cres[ind],
    subclass = ptsc,
    ATAC = atac[ind],
    H3K27ac = k27[ind]
  )
}) |>
  do.call(what = rbind, args = _)

scEnhs <- scEnhs[scEnhs$H3K27ac > 0, ]

data.table::fwrite(x = scEnhs,
  file = file.path(workd, "out",
    "ambpt_subclassDARChrACREs_sortBy-ATAC-FromFig4D.csv"),
  col.names = T, row.names = F, sep = ","
)

# * further update based on ds1000
allEnhDARCREs <- readRDS(
  file.path(projd, "06.ChromHMM", "out",
    "distal.Chr-A.DE.ATACPeak.rds")
)

raweRNA <- readRDS(
  file.path(projd, "06.ChromHMM", "out",
    "eRNA.logRPKM.filteredDAREhn.pbysc.rds")
)
mPbyS <- raweRNA

ordscs <- names(allEnhDARCREs) # already ordered
ordscs <- ordscs[!ordscs %in% c("044_OB_Dopa_Gaba")]
getOrdCRE <- function(topk = 1000, eRNAfd = 1.1) {
  partCRElist <- lapply(ordscs, \(sc) {
    reduceCREbyeRNA(
      sc = sc, enhDARCREs, fdPbyS,
      topk = topk, eRNAfd = eRNAfd
    )
  })
  names(partCRElist) <- ordscs
  partCRElist
}

reduceCREbyeRNA <- function(sc, CRElist,
                            fdPbyS, topk = 2000, eRNAfd = 1.1) {
  message(str_glue("reduce CREs by eRNA: {sc}."))
  m <- CRElist[[sc]]
  f <- fdPbyS[m, sc]
  index <- (f >= eRNAfd)
  if (sum(index) < 1) {
    message(str_glue("{sc} has no CREs under eRNA fd {eRNAfd}."))
    return(NULL)
  }
  mm <- m[index]
  return(mm[seq_len(min(topk, sum(index)))])
}

retrieveCRE <- function(sc, deATAC, log2fd = 0.2, n = 99999999L) {
  m <- deATAC[[sc]]
  index <- m[, 2] >= log2fd
  if (sum(index) < 1) {
    message(str_glue("{sc} has no CREs under log2fd {log2fd}."))
    index <- seq_len(nrow(m))
  }
  r <- m[order(m[index, 2], decreasing = T), 1]
  return(r[seq_len(min(length(r), n))])
}

enhDARCREs <- lapply(ordscs, retrieveCRE,
  deATAC = allEnhDARCREs,
  log2fd = 0.2, n = 999999
) |>
  setNames(object = _, nm = ordscs)
fdPbyS <- data.table::fread(
  file = file.path(projd, "06.ChromHMM", "out", "fd.csv"),
  header = F, sep = ",", data.table = F
) |>
  setNames(object = _, nm = colnames(mPbyS)) |>
  x => `rownames<-`(x, rownames(mPbyS))
partCREs <- getOrdCRE(topk = 1000, eRNAfd = 1.1)

scEnh2 <- scEnhs[scEnhs$subclass %in% ordscs, ]
scEnh3 <- lapply(ordscs, \(sc) {
  s <- scEnh2[scEnh2$subclass == sc, ]
  cres <- partCREs[[sc]][partCREs[[sc]] %in% s$cCRE]
  ind <- match(s$cCRE, cres) |>
    x => x[!is.na(x)]
  s[ind, ]
}) |>
  do.call(what = rbind, args = _)

data.table::fwrite(x = scEnh3,
  file = file.path(workd, "out",
    "ambpt_subclassDARChrACREs_sortBy-ATAC-FromFig4D.csv"),
  col.names = T, row.names = F, sep = ","
)
