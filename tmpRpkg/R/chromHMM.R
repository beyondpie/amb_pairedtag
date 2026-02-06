#' Read dense bed file from ChromHMM.
#' Ignore first row since it's a title for showing in track.
#'
#' @param densefnm chromHMM segmentation dense fnm
#' @return data.frame
#' @export
readChromHMMDenseBed <- function(densefnm) {
  r <- data.table::fread(densefnm,
    header = F, skip = 1, sep = "\t",
    data.table = F
  ) |> setNames(
    object = _,
    nm = c(
      "chr", "start", "end", "state",
      "score", "strand", "start2", "end2", "RGB"
    )
  )
  ## match GenomicRanges convention
  r$strand <- "*"
  return(r)
}


#' Read ChromHMM Binarized Data
#'
#' @param prefix path prefix
#' The binarized file is typically named as
#' prefix_[chr]_binary.txt
#' @param chrs vector of chromosomes for to read
#' @return list of data.frames with columns
#' id (start from 1) followed by different modules.
#' The names of the data.frames are chrs. If some chrs have no
#' binary files, will skip the corresponding chrs.
#' If all the chrs have no binary files, return NULL.
#' @export
loadChromHMMBinarizedData <- function(prefix, chrs) {
  lapply(chrs, \(chr) {
    fnm <- str_glue("{prefix}_{chr}_binary.txt")
    if (!file.exists(fnm)) {
      message(fnm, " does not exist, and skip it.")
      return(NULL)
    }
    r <- data.table::fread(
      file = fnm,
      header = T, skip = 1, data.table = F
    )
    r$id <- seq_len(nrow(r))
    rownames(r) <- paste0("id", r$id)
    return(r)
  }) |>
    setNames(object = _, nm = chrs) |>
    filterNULLfromList(l = _)
}

#' Map genomic bins to chromHMM states.
#' @param segment data.frame genomic states read from chromHMM
#' @param chrom2Bin data.frame from paritionChrom for one chromsome
#' @param ignoreStrand bool default TRUE sine chrom2Bin has no
#' strand information.
#' @return data.frame with the columns id (chromBin id) and state
#' defined in chromHMM. The rownames are "id1", "id2", ...
#' @export
mapBin2State.Chrom <- function(segment, chrom2Bin,
                               ignoreStrand = T) {
  chrom2BinGR <- loadGRfromBed(beds = chrom2Bin, header = T)
  segmentGR <- loadGRfromBed(beds = segment, header = T)
  hits <- GenomicRanges::findOverlaps(
    query = chrom2BinGR,
    subject = segmentGR,
    ignore.strand = ignoreStrand
  )
  r <- data.frame(
    id = chrom2Bin[S4Vectors::queryHits(hits), "id"],
    state = segment[S4Vectors::subjectHits(hits), "state"]
  )
  r <- r[!duplicated(r$id), ]
  rownames(r) <- paste0("id", r$id)
  return(r)
}

#' Get the emission counts for one chromosome.
#' @param bin2PromoterDistal data.frame from the function
#' mapBin2PromoterDistal
#' @param bin2State data.frame from the function mapBin2State.Chrom
#' @param bin2BinarizeData data.frame from the function
#' loadChromHMMBinarizedData and focus on one chromsom.
#' @return data.frame with id, state, relaTSS, pdstate,
#' different modularites and gsize.
#' Add all the counts for those modularities.
#' gsize is the total counts for each of the pdstate.
#' @export
getEmissionCount.Chrom <- function(bin2PromoterDistal,
                                   bin2State,
                                   bin2BinarizeData) {
  # relaTSS ("P" or "Q"), State (integer defined from ChromHMM),
  # different modularities
  r <- bin2PromoterDistal |>
    dplyr::inner_join(x = _, y = bin2State, by = c("id")) |>
    dplyr::inner_join(x = _, y = bin2BinarizeData, by = c("id"))
  r$pdstate <- with(r, paste(state, relaTSS, sep = "_"))
  rr <- r |>
    dplyr::group_by(pdstate) |>
    summarise(
      .data = _,
      across(
        -any_of(c("id", "state", "relaTSS", "pdstate")),
        \(x) sum(x, na.rm = TRUE)
      ),
      gsize = n()
    ) |>
    as.data.frame(x = _) |>
    x => `rownames<-`(x, x$pdstate)
  return(rr)
}

#' Use the list of EmissionCounts for different chromosomes
#' to estimate the emission probabilisties from the learned data.
#'
#' @param listOfEmissionCount list of emission counts from the function
#' of getEmissionCount.Chrom
#' @param mods vector of modularities.
#' @param states vector of states
#' @return data.frame with columns of
#' raw counts per mod, totalcount for states, and the avgMat
#' per mod.
#' @export
calAvgEmissionMat <- function(listOfEmissionCount,
                              mods = c(
                                "ATAC",
                                "H3K27ac", "H3K27me3", "H3K4me1",
                                "H3K9me3"
                              ),
                              states = ordStates) {
  # fix states not exists
  l <- lapply(listOfEmissionCount, \(t) {
    t <- t[, c("pdstate", mods, "gsize")]
    if (!all(states %in% rownames(t))) {
      missStates <- setdiff(states, rownames(t))
      t_left <- matrix(0,
        nrow = length(missStates), ncol = ncol(t) - 1,
        dimnames = list(missStates, c(mods, "gsize"))
      ) |> as.data.frame()
      t_left$pdstate <- missStates
      t <- rbind(t, t_left)
    }
    rownames(t) <- t$pdstate
    t
  })
  # align all the matrix
  alignCountsOfMods <- lapply(l, \(t) {
    t[states, mods] |> as.matrix()
  }) |>
    Reduce(`+`, x = _)

  alignTotalCounts <- lapply(l, \(t) {
    t[states, "gsize"]
  }) |> Reduce(`+`, x = _)

  # each row is divided by the same n_total in the total counts.
  r <- alignCountsOfMods / (alignTotalCounts + 1e-10)
  r <- as.data.frame(r) |>
    x => `rownames<-`(x, states)

  data.frame(
    countMat = alignCountsOfMods,
    totalCount = alignTotalCounts,
    avgMat = r
  )
}
