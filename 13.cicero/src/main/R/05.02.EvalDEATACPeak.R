library(tidyverse)
detach("package:tmpRpkg", unload = T)
library(tmpRpkg)

# * meta
stateNameOrd <- tmpRpkg:::stateNameOrd
projd <- here::here()
ATACded <- file.path(projd, "12.DE", "out", "ATAC_sa2LRT_DE")

ptscMeta <- tmpRpkg::loadptscMeta() |>
  x => ptscMeta[x$ATAC > 0, ]
scsChromHMM <- ptscMeta$PairedTagName

scsDE <- list.files(ATACded, include.dirs = F, no.. = T) |>
  lapply(X = _, \(f) gsub("_ATAC_sa2DE.tsv", "", f)) |>
  unlist()

## homer info
mm10chromSizefnm <- file.path(
  projd,
  "meta", "mm10.chrom.sizes.lite"
)
chromSizes <- data.table::fread(
  mm10chromSizefnm,
  header = F, data.table = F
) |>
  setNames(object = _, nm = c("chr", "size")) |>
  x => `rownames<-`(x, x$chr)

## use unified promoter: TSS +- 1kb
promoter <- tmpRpkg::loadPromoter()

## chromHMM annot
annotd <- file.path(projd, "06.ChromHMM", "out", "CREAnnot250429")

## output
outd <- file.path(projd, "13.cicero", "out")
figd <- file.path(projd, "13.cicero", "figure")

# * functions
loadDE <- function(fnm, q = 0.05) {
  data.table::fread(file = fnm, sep = "\t", header = T, data.table = F) |>
    x => `rownames<-`(x, x[, 1]) |>
    x => x[x[, 4] <= q, ]
}

# * main
# now all the sc from chromHMM has ATAC DE result.
# scs <- intersect(scsDE, scsChromHMM)
scs <- scsChromHMM
qval <- 0.05
state <- "Chr-A"
# No Chr-P observed in oligo, ast.
activeCREscs <- lapply(scs, \(sc) {
  DEs <- loadDE(
    file.path(ATACded, str_glue("{sc}_ATAC_sa2DE.tsv")),
    q = qval
  )
  message(str_glue("{sc} has {nrow(DEs)} DE CREs."))

  rawAnnot <- readSummitAnnot(
    scfnm = file.path(annotd, str_glue("{sc}.CRE.annotBySummit.bed")),
    promoter,
    sc = sc)
  unstate <- length(unique(rawAnnot$state))
  if (unstate > 7) {
    message(str_glue("{sc} has {unstate} kinds of distal CREs."))
  }

  annotCREs <- rawAnnot |>
    x => x[x$state == state, ] |>
    x => x[x$pd == "distal", ]

  message(str_glue("{sc} has {nrow(annotCREs)} distal {state} CREs."))

  r <- intersect(annotCREs[, 1], DEs[, 1])
  message(str_glue("{sc} has {length(r)} specfic CREs of distal and {state}."))
  DEs[r, ]
}) |>
  setNames(object = _, nm = scs)

saveRDS(
  object = activeCREscs,
  file = file.path(outd, str_glue("distal.{state}.DE.ATACPeak.rds"))
)

allDEDistalCREs <- lapply(scs, \(sc) {
  DEs <- loadDE(
    file.path(ATACded, str_glue("{sc}_ATAC_sa2DE.tsv")),
    q = qval
  )
  message(str_glue("{sc} has {nrow(DEs)} DE CREs."))

  rawAnnot <- readSummitAnnot(
    scfnm = file.path(annotd, str_glue("{sc}.CRE.annotBySummit.bed")),
    promoter,
    sc = sc)
  
  unstate <- length(unique(rawAnnot$state))
  if (unstate > 7) {
    message(str_glue("{sc} has {unstate} kinds of CREs."))
  }

  annotCREs <- rawAnnot |>
    x => x[x$pd == "distal", ]

  message(str_glue("{sc} has {nrow(annotCREs)} distal CREs."))

  r <- intersect(annotCREs[, 1], DEs[, 1])
  message(str_glue("{sc} has {length(r)} specfic distal CREs."))
  DEs[r, ]
}) |>
  setNames(object = _, nm = scs)

saveRDS(
  object = allDEDistalCREs,
  file = file.path(outd, str_glue("all.distal.DE.ATACPeak.rds"))
)
