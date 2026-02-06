#' Transform Bed to GenomicRanges
#' Bed is from a bedFile or a data.frame
#'
#' @param bedFile str or NULL
#' @param beds data.frame or NULL
#' for data.frame, at least 3 columns
#' (first 3 columns are chrom, start end)
#' @param header boolean, default FALSE
#' @param sep str, default tab
#' @param colnms List[str], default is
#' c("chr", "start", "end", "name")
#' @return GenomicRanges
#' @export
loadGRfromBed <- function(bedFile = NULL,
                          beds = NULL,
                          header = FALSE,
                          sep = "\t",
                          colnms = c("chr", "start", "end", "name")) {
  if (is.null(beds) & is.null(bedFile)) {
    stop("At lease one of bedFile or beds should be NOT NULL.")
  }
  if (is.null(beds)) {
    beds <- data.table::fread(
      file = bedFile,
      sep = sep, header = header
    )
  }
  if (!header) {
    if (ncol(beds) < 4) {
      ## for beds without names
      colnms <- colnms[1:3]
    }
    colnames(beds) <- colnms
  } else {
    colnms <- colnames(beds)
  }

  gr <- GenomicRanges::GRanges(
    seqnames = beds[[1]],
    ranges = IRanges::IRanges(
      start = beds[[2]],
      end = beds[[3]]
    )
  )
  if ("strand" %in% colnms) {
    rawStrand <- beds$strand
    idx <- !(rawStrand %in% c("+", "-", "*"))
    if (sum(idx) > 0) {
      message("Characters out of +, -, * are found in strand column.")
      message("Treat them as *.")
      rawStrand[idx] <- "*"
    }
    GenomicRanges::strand(gr) <- rawStrand
  }
  unWantedNames <- c(
    "seqnames", "ranges", "strand",
    "seqlevels", "seqlengths", "isCircular", "start",
    "end", "width", "element"
  )

  if (ncol(beds) >= 4) {
    for (j in seq(from = 4, to = ncol(beds))) {
      nm <- colnms[[j]]
      if (nm != "strand") {
        if (nm %in% unWantedNames) {
          message(nm, "is not recommended in GenomicRanges.")
          nm <- paste0(nm, ".1")
          message("Change it into ", nm)
        }
        S4Vectors::mcols(gr)[, nm] <- beds[[j]]
      }
    }
  }
  return(gr)
}

#' Transform GenomicRanges Object to chr:start=end format
#' @param gr GenomicRanges Object
#' @return Array of characters
#' @export
strGRSeq <- function(gr) {
  chrs.tmp <- GenomicRanges::seqnames(gr)
  chrs.tmp <- as.factor(chrs.tmp)
  chrs <- levels(chrs.tmp)[chrs.tmp]
  starts <- GenomicRanges::start(gr)
  ends <- GenomicRanges::end(gr)
  r <- paste(chrs, paste(starts, ends, sep = "-"), sep = ":")
  return(r)
}

#' Transform GenomicRanges Object to
#' three-column bed (simple bed)
#' @param gr GenomicRanges object
#' @return data.frame three columns
#' @export
transformGR2SimpleBed <- function(gr) {
  chrs.tmp <- GenomicRanges::seqnames(gr)
  chrs.tmp <- as.factor(chrs.tmp)
  r <- data.frame(
    chr = levels(chrs.tmp)[chrs.tmp],
    startFrom = GenomicRanges::start(gr),
    endTo = GenomicRanges::end(gr)
  )
  return(r)
}

#' Modification of GenomicRanges::findOverlaps by adding
#' overlap region, width.
#' @export
findOverlaps2 <- function(queryGR,
                          subjectGR) {
  hits <- GenomicRanges::findOverlaps(
    query = queryGR, subject = subjectGR
  )
  queryIndex <- S4Vectors::queryHits(hits)
  subjectIndex <- S4Vectors::subjectHits(hits)
  ## get intersect regions
  s <- GenomicRanges::pintersect(
    x = queryGR[queryIndex],
    y = subjectGR[subjectIndex]
  )
  s$ovlpWidth <- GenomicRanges::width(s)
  ## add subject and query seqs
  s$queryGR <- queryGR[queryIndex]
  s$subjectGR <- subjectGR[subjectIndex]
  ## restore hits info
  s$hit <- NULL # remove original TRUE label
  s$ovlp <- hits
  return(s)
}

#' @export
filterByGR <- function(A, B) {
  .Deprecated("Use FindOvlpRegionInB instead.")
  FindOvlpRegionInB(A, B)
}

#' Find the genomic regions in A that overlap with B
#'
#' @export
FindOvlpRegionInB <- function(A, B) {
  hits <- GenomicRanges::findOverlaps(query = A, subject = B)
  return(A[unique(S4Vectors::queryHits(hits))])
}

#' tranformRegion2GR
#' @param region vector of "chr1:1-2" format
#' @return data.frame with the column names
#' c(chr, startFrom, endTo, name)
#' @export
transformRegion2GR <- function(region) {
  strsplit(region, split = ":|-", perl = T) |>
    do.call(rbind, args = _) |>
    as.data.frame() |>
    setNames(object = _, nm = c("chr", "startFrom", "endTo")) |>
    mutate(
      startFrom = as.integer(startFrom),
      endTo = as.integer(endTo),
      name = region
    ) |> loadGRfromBed(beds = _)
}


#' transformGenomicRange2DataFrame
#' @export
transformGenomicRange2DataFrame <- function(
    gr, zeroBased = TRUE, default_strand = ".") {
  outdf <- as.data.frame(gr)
  for (i in colnames(outdf)) {
    if (is.factor(outdf[, i])) {
      # message("Transform factor to characters: ", i)
      outdf[, i] <- levels(outdf[, i])[outdf[, i]]
    }
  }
  outdf$strand <- gsub("\\*", default_strand, outdf$strand)
  return(outdf)
}

#' From GTF (loaded using rtracklayer::import from gtf / gff[3])
#' GenomicRanges, get the gene GenomeRanges based on its type
#' @param gtf GenomicRanges
#' @param mcolnm vector of characters metadata for result
#' @return GenomicRanges
#' @export
getGeneGR <- function(gtf,
                      mcolnms = c("gene_id", "type", "gene_type", "gene_name")) {
  types <- S4Vectors::mcols(gtf)$type
  types <- levels(types)[types]
  index.gene <- types %in% "gene"
  geneGTF <- gtf[index.gene, mcolnms]
  return(geneGTF)
}

#' @export
aroundGeneTSS <- function(gtf,
                          up = 1000,
                          down = 1000,
                          mcolnms = c("gene_id", "type", "gene_type", "gene_name"),
                          ...) {
  types <- S4Vectors::mcols(gtf)$type
  types <- levels(types)[types]
  index.gene <- types %in% "gene"
  geneGTF <- gtf[index.gene, mcolnms]
  r <- GenomicRanges::promoters(
    geneGTF,
    upstream = up, downstream = down, ...
  )
  return(r)
}

#' @export
getGeneTSS <- function(gtf,
                       mcolnms = c("gene_id", "type", "gene_type", "gene_name"),
                       ...) {
  r <- aroundGeneTSS(gtf = gtf, up = 0, down = 1, mcolnms = mcolnms, ...)
  return(r)
}

#' @export
resizeGeneTSS <- function(gtf,
                          window.size = 1000,
                          mcolnms = c("gene_id", "type", "gene_type", "gene_name"),
                          fix = "center",
                          ...) {
  tss <- getGeneTSS(gtf = gtf, mcolnms = mcolnms)
  r <- GenomicRanges::resize(x = tss, width = window.size, fix = fix, ...)
  return(r)
}

#' Resize the GenomicRanges with ext_size and ignore the strand info.
#'
#' @param gr GenomicRanges
#' @param fix str, can be c("center", "start", "end", "twoside")
#' @param extsize int, default 2500
#' @return GenomicRanges
#' @export
resizeIgnoreStrand <- function(gr, fix = "center", extsize = 2500) {
  if (fix == "twoside") {
    wids <- GenomicRanges::width(gr)
    r1 <- GenomicRanges::resize(gr,
      width = wids + extsize - 1, fix = "start",
      ignore.strand = TRUE
    )
    r2 <- GenomicRanges::resize(gr,
      width = extsize + wids - 1, fix = "end",
      ignore.strand = TRUE
    )
    r <- gr
    GenomicRanges::start(r) <- GenomicRanges::start(r2)
    GenomicRanges::end(r) <- GenomicRanges::end(r1)
    return(r)
  }
  return(GenomicRanges::resize(gr,
    width = extsize * 2,
    fix = fix, ignore.strand = TRUE
  ))
}


#' Given one chrom size, partition it into bins.
#'
#' @param chrName one chromosome name
#' @param chromSize integer size of the chrom
#' @param binSize integer
#' @param excludeLastPartialBin default TRUE
#' ChromHMM will exclude the last bin if it's small than binSize
#' @return data.frame of chr, start, end, bin id (from 1)
#' chrom labels following bed format, i.e.,
#' start from 0, and left closed and right open
#' for example [0, 1), 1st position.
#' @export
partitionChrom <- function(chrName, chromSize,
                           binSize = 200,
                           excludeLastPartialBin = T) {
  binStarts <- seq(0, chromSize - 1, by = binSize)
  binEnds <- pmin(binStarts + binSize, chromSize)
  r <- data.frame(
    chr = chrName,
    start = binStarts,
    end = binEnds,
    id = seq_along(binStarts)
  )
  rownames(r) <- paste0("id", r$id)
  # check last bin is full size or not
  lastRow <- nrow(r)
  if ((binEnds[lastRow] - binStarts[lastRow]) < binSize) {
    if (excludeLastPartialBin) {
      r <- r[-lastRow, ]
    }
  }
  return(r)
}


#' Map Genomic bins to promoer or distals based on homer TSS.
#'
#' @param chrom2Bin data.frame got from paritionChrom
#' Here chrom2Bin is for one chromosome.
#' @param homerTSS data.frame got from readHomerTSS
#' @param upSize integer default 1000 up to TSS
#' @param downSize integer default 1000 down to TSS
#' @param ignoreStrand bool default TRUE
#' chrom2Bin has no strand information.
#' @return data.frame with the columns of
#' id, chrom2Bin's id; relaTSS, "D" for distals, "P" for
#' promoters.
#' The rownames are "id1", "id2"...
#' @export
mapBin2PromoterDistal.Chrom <- function(chrom2Bin, homerTSS,
                                        upSize = 1000,
                                        downSize = 1000,
                                        ignoreStrand = T) {
  chrom2BinGR <- loadGRfromBed(beds = chrom2Bin, header = T)
  homerTSSGR <- loadGRfromBed(beds = homerTSS, header = T)
  homerPromoter <- GenomicRanges::promoters(
    x = homerTSSGR, upstream = upSize, downstream = downSize
  )
  hits <- GenomicRanges::findOverlaps(
    query = chrom2BinGR, subject = homerPromoter,
    ignore.strand = ignoreStrand
  )
  r <- data.frame(
    id = chrom2Bin$id,
    relaTSS = "D",
    row.names = paste0("id", chrom2Bin$id)
  )
  ## rows located in the query Hits are proximals.
  r$relaTSS[
    seq_len(nrow(r)) %in% (S4Vectors::queryHits(hits))
  ] <- "P"
  return(r)
}

#' Extend upstream and downsteam of Genomic ranges
#' 
#' @param g genomic ranges
#' @param up upstream extSize, positive integer
#' @param down downstream extSize, positive integer
#' @return a new Genomic ranges with extSizes.
#' @export
extendGR <- function(g, up = 500, down = 500) {
  r <- g
  GenomicRanges::start(r) <- GenomicRanges::start(g) - up
  GenomicRanges::end(r) <- GenomicRanges::end(g) + down
  return(r)
}
