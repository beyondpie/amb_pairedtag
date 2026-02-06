#' @export
transformAllenLabel2str <- function(allenLabels) {
  vapply(allenLabels, \(a) {
    gsub(" |-|/", "_", a)
  }, "a")
}

#' @export
sortAllenLabelById <- function(allenLabels) {
  ids <- vapply(allenLabels, \(a) {
    as.integer(stringr::str_split_1(a, "_")[1])
  }, 1)
  allenLabels[order(ids)]
}

#' @export
load_sp2sc <- function() {
  allenClMeta <- loadAllenClMeta()
  sp2sc <- allenClMeta[,
    c("supertype_id_label", "subclass_id_label")] |>
    unique() |>
    setNames(object = _, nm = c("sp", "sc"))
  sp2sc$sp <- transformAllenLabel2str(sp2sc$sp)
  sp2sc$sc <- transformAllenLabel2str(sp2sc$sc)
  rownames(sp2sc) <- sp2sc$sp
  return(sp2sc)
}

#' Map subclass in snATAC name habit to pt subclass id label habit
#' @return data.frame with sc_atac and sc as two columns
#' @export
load_sc2scil <- function() {
  allenClMeta <- loadAllenClMeta()
  sc2scil <- allenClMeta[, c("subclass_label", "subclass_id_label")] |>
    unique() |>
    setNames(object = _, nm = c("sc_atac", "sc"))
  sc2scil$sc_atac <- stringr::str_replace_all(
    sc2scil$sc_atac, " +", "_") |>
    x => stringr::str_replace_all(x, "/", "-")
  sc2scil$sc <- transformAllenLabel2str(sc2scil$sc)
  rownames(sc2scil) <- sc2scil$sc_atac
  return(sc2scil)
}

#' @export
load_snATACPMpbyc <- function() {
  r <- data.table::fread(snATACPMpbyc, header = TRUE,
    data.table = FALSE)
  rownames(r) <- r$V1
  r$V1 <- NULL
  old_names <- colnames(r)
  sc2scil <- load_sc2scil()
  colnames(r) <- sc2scil[old_names, "sc"]
  return(r)
}

#' Perform CPM for Histone raw feature counts.
#' @param cnt_pbyc matrix, feature by group
#' @param scale.factor default 1e6
#' @return CPM, feature by group
#' @export
getCPMOfHistones <- function(cnt_pbyc, scale.factor = 1e6) {
  t(t(cnt_pbyc * scale.factor) / colSums(cnt_pbyc))
}

#' Get CPM normalized Paired-Tag gene expressions.
#'
#' BUG:
#' - The original data is data.frame, which has to been converted into
#'   matrix. Otherwise, will have weired behaviour on CPM normalization.
#'   (i.e., matrix / vector and data.frame / vector are different?)
#' - "|>" piple operator leads to error during noramlization.
#'   + Error in (function(x) (x * 10^6)) ((function(x) x[atac2ptsc$pt, ])
#'   (as.matrix(readRDS(ptPseudoRNAcountfnm)))) %*%
#'   :non-conformable arguments
#' 
#' @param fnm \link{ptPseudoRNAcountfnm}
#' R's rds file, subclass by gene expression.
#' @param scaleFactor used to scaling, default is 1e6 (million)
#' @return matrix, group by feature
#' @export
getCPMOfptRNA <- function(fnm, scaleFactor = 1e6) {
  m <- readRDS(fnm) |>
    as.matrix()
  (m * scaleFactor) / rowSums(m)
}

#' Map Histone raw feature counts to subclass-level
#' Non-neuronal groups are supertypes, so we need a map.
#' @export
mapfc2fcsc <- function(fc, sp2scMap) {
  mapsp2scMat <- function(spmat, sp2scMap) {
    sps <- colnames(spmat)
    scs <- sp2scMap[sps, "sc"]
    uscs <- unique(scs)
    r <- vapply(uscs, \(sc) {
      ids <- scs %in% sc
      if (sum(ids) < 2) {
        spmat[, ids]
      } else {
        rowSums(spmat[, scs %in% sc])
      }
    }, spmat[, 1])
    colnames(r) <- uscs
    return(r)
  }
  fcsc1 <- fc[, colnames(fc) %in% sp2scMap$sc]
  fcsp <- fc[, colnames(fc) %in% sp2scMap$sp]
  fcsc2 <- mapsp2scMat(fcsp, sp2scMap)
  cbind(fcsc1, fcsc2)
}

#' map subclass name to snATAC subclass format
#' @param sc vector of string
#' @return vector of string
#' @export
map2snATACscstr <- function(sc) {
  str_replace_all(sc, "^\\d+ ", "") |>
    x => str_replace_all(x, " +",  "_") |>
    x => str_replace_all(x, "/", "-") 
}
