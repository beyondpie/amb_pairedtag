library(GenomicRanges)
library(BSgenome)
gr <- GRanges(Rle(c("chr2", "chr1", "chr3", "chr1"), 4:1),
                  IRanges(1:10, width=5),
                  seqinfo=Seqinfo(c("chr1", "chr2", "chr3"), c(100, 50, 20)),
                  score = rnorm(n = 10, sd = 10, mean = 5)
) |>
  sortSeqlevels() |>
  sort()
r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
o <- findOverlaps(gr,r, ignore.strand = TRUE)
mcols(gr)$cluster <- subjectHits(o)
gr <- gr[order(mcols(gr)[, "score"], decreasing = TRUE),]
grSelect <- gr[!duplicated(mcols(gr)$cluster),] |>
  sortSeqlevels() |>
  sort()
grConverge <- subsetByOverlaps(
  gr,
  grSelect, 
  invert=TRUE, 
  ignore.strand = TRUE) #blacklist selected gr

genome <- getBSgenome("mm10")
nucFreq <- BSgenome::alphabetFrequency(getSeq(genome, gr))

# * help check scala's implementation of overlap
library(rtracklayer)
library(tidyverse)
library(GenomicRanges)
library(tmpRpkg)

projd <- here::here()
workd <- file.path(projd, "13.cicero")
ptRNAbwd <- file.path(projd, "data", "ptRNAbam", "bigwig")
atacPeakd <- file.path(projd, "data", "chromHMM", "subclass_peak")

loadptRNAbw <- function(sc) {
  rtracklayer::import(con = file.path(
    ptRNAbwd, str_glue("{sc}.RPKM.bw")
  ))
}

loadscATAC <- function(sc, chrs = c("chr1")) {
  data.table::fread(
    file = file.path(atacPeakd, str_glue("{sc}.ATAC.peak.bed")),
    header = F, sep = "\t", data.table = F
  ) |>
    setNames(object = _, nm = c("chr", "start", "end", "name")) |>
    x => x[x$chr %in% chrs, ] |>
    loadGRfromBed(beds = _, header = F)
}

sc <- "334_Microglia_NN"
chr <- "chr1"
bw <- loadptRNAbw(sc)
CREs <- loadscATAC(sc)

# * test ovlp
ovlps <- findOverlaps(query = CREs, subject = bw)
qids <- ovlps@from-1
sids <- ovlps@to-1
r <- data.frame(
  q = qids,
  s = sids
) |> group_by(q) |>
  summarise(
    sleft = min(s),
    sright = max(s) + 1)

# * test scores
chr <- "chr1"
geteRNAlogRPKM <- function(sc, CREs) {
  ptRNAbw <- loadptRNAbw(sc)
  message(str_glue("get {sc}'s eRNA RPKM."))
  pteRNAall <- getbwSignal(bwGR = ptRNAbw, peakGR = CREs)
  data.frame(
    CRE = names(pteRNAall),
    pteRNA = pteRNAall
  )
}

scores <- geteRNAlogRPKM(sc, CREs) |>
  x => x[grepl(chr, x$CRE), ]
