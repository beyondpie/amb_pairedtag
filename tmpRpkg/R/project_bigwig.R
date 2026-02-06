pteRNAbwd <- file.path(projd, "data", "ptRNAbam", "bigwig")
ptDNAbwd <- file.path(projd, "data", "ptDNAbam", "bigwig")

#' Load one histone's bigwig as Genomic ranges
#' @param sc PairedTag subclass name
#' @export
loadHistoneBigWig <- function(sc, h) {
  suffix <- if (h == "H3K9me3") {
    "e100.bs100.sm1000"
  } else {
    "e100.bs100.sm300"
  }
  fnm <- file.path(ptDNAbwd, stringr::str_glue(
    "{sc}.{h}.{suffix}.bw"
  ))
  loadBigWigGR(fnm)
}

