# * global variables

snATACd <- file.path(projd, "data", "snATAC")
snATACSrtPeakbedfnm <- file.path(snATACd, "wmb.snATAC.peak.srt.bed")
snATACPeakInPt <- file.path(snATACd, "snATAC.peaks.included.pairedtag.csv")
snATACPMpbyc <- file.path(snATACd, "cpm_peakBysubclass.csv")
snATACountpbyc <- file.path(snATACd, "count_peakBysubclass.csv")
snATACProjd <- "/projects/ps-renlab2/szu/projects/CEMBA2/"
snATACsa2DEd <- file.path(projd, "12.DE", "out", "ATAC_sa2LRT_DE")
snATAChrADEfnm <- file.path(projd, "13.cicero", "out",
  "distal.Chr-A.DE.ATACPeak.rds")
snATACDEfnm <- file.path(projd, "13.cicero", "out",
  "all.distal.DE.ATACPeak.rds")
snATACppecd <- file.path(projd, "13.cicero", "out", "ATACppec")

#' Get ATAC Peaks' names
#' @return vector of string
#' @export
loadsrtCREname <- function() {
  data.table::fread(
    file = snATACSrtPeakbedfnm, header = FALSE, data.table = FALSE
  )$V4
}

#' Load All one million CREs from ATAC-seq data as Genomic Ranges
#' @export
loadAllCREAsGR <- function() {
  loadGRfromBed(bedFile = snATACSrtPeakbedfnm)
}

#' Get ATAC Peaks' located in PairedTag region
#' @return vector of string
#' @export
loadsrtCREnameInPairedTag <- function() {
  p <- loadsrtCREname()
  data.table::fread(
    # snATACPeakInPt 
    file = snATACPeakInPt, header = FALSE, data.table = FALSE
  )$V1 |>
    x => p[p %in% x]
}
#' Get chromHMM annotations for ATAC-seq peaks.
#'
#' @param scfnm state-mapped ATAC-summit file.
#' Under [projd]/06.ChromHMM/out/CREAnnot250429
#' the filename is like [subclass].CRE.annotBySummit.bed
#' @param promoter GRobject  see \link{getHomerPromoterGR} as example
#' @param sc subclass name
#' @return data.frame with columns CRE, state,
#' pd (short for proximal or distal)
#' @export
readSummitAnnot <- function(scfnm, promoter, sc = NULL) {
  r <- data.table::fread(
    file = scfnm,
    sep = "\t",
    header = F,
    data.table = F
  ) |>
    setNames(object = _,
      nm = c("chrom", "startFrom", "endTo", "state"))
  pGR <- loadGRfromBed(beds = r)
  hits <- GenomicRanges::findOverlaps(query = pGR, subject = promoter)
  proxIndex <- S4Vectors::queryHits(hits) |>
    sort() |>
    unique()
  pd <- rep("distal", length(pGR))
  pd[proxIndex] <- "proximal"
  w <- data.frame(
    CRE = with(r, paste(chrom,
      paste(startFrom, endTo, sep="-"), sep = ":")),
    state = r$state,
    pd = pd,
    sc = sc
  ) |> x => `rownames<-`(x, x$CRE)
  return(w)
}


#' Load ATAC pseudo-bulk CPM / Count as matrix.
#' @param f \link{snATACPMpbyc} or \link{snATACountpbyc}
#' @param pd pandas module from reticulate::import
#' if NULL, will set up it here.
#' @param env conda env used for pd, "sa2" by default
#' @param conda conda bin path
#' @return matrix with both rownames and colnames.
#' @export
loadATACPseudoBulkMat <- function(f, pd = NULL,
                                  env = "sa2",
                                  conda =
                                    "/home/szu/mambaforge/bin/conda") {
  if (is.null(pd)) {
    message(str_glue(
      "set up python pandas env: conda-{conda} env-{env}."
    ))
    tmpRpkg::setCondaEnv(envnm = env, conda = conda)
    pd <- reticulate::import(module = "pandas", convert = F)
  }
  pd$read_csv(f, index_col = 0L, sep = ",") |>
    reticulate::py_to_r(x = _) |>
    as.matrix()
}

#' Filter deATAC based on log2 fold change
#'
#' @param sc PairedTag-format subclass
#' @param deATAC \link{snATAChrADEfnm}
#' The CREs has already be filtered by q-value 0.05, though
#' we should do it here for raw DE files.
#' @param log2fd 0.2 by default
#' @param n number of top ranked CREs
#'
#' @return vector of CREs (chr1:0-1 format)
#' CREs will be sorted based on log2fd
#'
#' @export
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
