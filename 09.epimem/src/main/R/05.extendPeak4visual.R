Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

# * meta
extsize <- 5000

# * set arguments
args <- commandArgs(trailingOnly = TRUE)
infnm <- args[1]
outfnm <- args[2]

extractBedFromGR <- function(gr) {
  chrs.tmp <- GenomicRanges::seqnames(gr)
  chrs.tmp <- as.factor(chrs.tmp)
  chrs <- levels(chrs.tmp)[chrs.tmp]
  starts <- GenomicRanges::start(gr)
  ends <- GenomicRanges::end(gr)
  r <- data.frame(
    chrom = chrs,
    start = starts,
    end = ends
  )
  return(r)
}

# * main
gr <- tmpRpkg::loadGRfromBed(bedFile = infnm, header = TRUE, sep = "\t")

r <- tmpRpkg::resizeIgnoreStrand(gr, fix = "twoside", extsize = extsize) |>
  extractBedFromGR()

data.table::fwrite(x = r, file = outfnm, col.names = FALSE, sep="\t")
