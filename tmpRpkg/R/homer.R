#' Load the known motifs learned from Homer2.
#' The file end with knownResults.txt will be handled here.
#' @param resultd where the homer2 results are saved
#' @param rmdupmotif default TRUE
#' @return data.frame
#' @export
loadHomer2KnownResults <- function(resultd, rmdupmotif) {
  # read results
  a <- data.table::fread(
    file = file.path(resultd, "knownResults.txt"),
    sep = "\t", header = T, data.table = F
  )
  # extract motif name from first column
  motifnms <- vapply(a[, 1], function(i) {
    unlist(strsplit(i, split = "/", fixed = TRUE))[1]
  }, FUN.VALUE = "CTCF")
  r <- data.frame(
    motif = motifnms,
    p = a[, 3],
    logp = a[, 4],
    qBenj = a[, 5],
    cnt = a[, 6],
    prcnt = a[, 7],
    cnt_bg = a[, 8],
    prcn_tbg = a[, 9]
  )
  ## temp remove duplicated motif
  if (rmdupmotif) {
    r <- r[!duplicated(r$motif), ]
  }
  return(r)
}

#' Summarize motif results across multiple groups
#' into a unified matrix.
#'
#' We assume that, for each motif result, it's a data.frame
#' following the formats defined in loadHomer2KnownResults function.
#'
#' @param motif_list list of homer results for different groups
#' @param field which column for the final matrix
#' @param tf function to transform the column
#' for example, we can use log10 or -log10 for negative
#' log10 of pvalue.
#' @return matrix groups as rowname, motifs as colnames
#' @export
getMotifMat <- function(motif_list,
                        field = "p",
                        tf = identity) {
  all.motif.nms <- unique(unlist(lapply(motif_list, \(i) i$motif)))
  # make sure all the motifs results have the same order of motifs
  a <- lapply(motif_list, \(i) {
    rownames(i) <- i$motif
    i[all.motif.nms, ]
  })
  r <- lapply(a, \(i) {
    tf(as.numeric(i[, field]))
  }) |> do.call(what = rbind, args = _)
  rownames(r) <- names(motif_list)
  colnames(r) <- all.motif.nms
  return(r)
}


#' Read Homer(i.e, Homer2) genome TSS file into data.frame.
#'
#' 1. Homer uses 2000 as buffersize to extend the raw TSS, which we
#'   will get the raw single base pair annotation.
#' 2. Homer use 1 as "-" strand and 0 as "+".
#' 3. The first column is RefSeq RNA id, so multiple isoforms shared
#'    the same TSS will also be listed.
#' Check "bin/parseGTF.pl" under the installed homer directory.
#' Check https://genome.ucsc.edu/FAQ/FAQformat.html for bed format.
#'
#' ChromHMM also provides the coords. Under instealled dir,
#' "COORDS/mm10/RefSeqTSS.mm10.bed.gz" as one example.
#' It only records the single-base, unique (no RefSeq id),
#' and NOT-ordered TSS. Both NR (non-coding) and NM (coding)
#' RNAs are recorded.
#'
#' @param homertssf full path of homer tss file.
#' For example, if the genome is mm10, it should be under
#' installed homer directory, "data/genomes/mm10/mm10.tss"
#' @return data.frame
#' @export
readHomerTSS <- function(homertssf, bufferSize = 2000) {
  r <- data.table::fread(
    file = homertssf,
    header = F, sep = "\t",
    data.table = F
  ) |>
    x => setNames(
      object = x,
      nm = c("RefSeqRNAid", "chr", "start", "end", "strand")
    ) |>
    dplyr::mutate(
      .data = _,
      start = as.integer(start - 1 + bufferSize),
      end = as.integer(end - bufferSize),
    )
  pos_strand <- r$strand == 0
  neg_strand <- r$strand == 1
  r$strand[pos_strand] <- "+"
  r$strand[neg_strand] <- "-"
  r <- r[, c("chr", "start", "end", "RefSeqRNAid", "strand")]
  return(r)
}

#' Get Promoter from Homer TSS
#' 
#' @param homertssf full path of homer TSS file
#' @param bufferSize 2000 See details in \link{readHomerTSS}
#' @param upstream 1000
#' @param downstream 500
#' @return GenomicRanges object
#' @export
getHomerPromoterGR <- function(homertssf,
                              bufferSize = 2000, upstream = 1000,
                              downstream = 500) {
  homerTSS <- readHomerTSS(homertssf, bufferSize) |>
    x => loadGRfromBed(beds = x, colnms = colnames(x))
  homerPromoter <- GenomicRanges::promoters(
    x = homerTSS, upstream = upstream, downstream = downstream)
  return(homerPromoter)
}


