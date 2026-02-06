# * meta
root <- "/projects/ps-renlab2/szu"
projd <- file.path(root, "projects/amb_pairedtag")
histones <- c("H3K27ac", "H3K4me1", "H3K27me3", "H3K9me3")
ptscMetafnm <- file.path(projd, "meta",
  "PairedTagSubclassMetaFromCellMetaChromHMM.csv")
hiscntd <- file.path(projd, "data", "ptHistoneCounts", "ATACPeak")
hiscntdr <- file.path(projd, "data", "ptHistoneCounts", "ATACPeak_rds")
ptPseudoRNAcountfnm <- file.path(projd, "data", "pairedtag_ann",
  "pt.RNAseq.pseudobulk.count.sc.allgenesymbol.rds")

homerd <- file.path(root, "softwares", "homer")
homermm10Annotd <- file.path(homerd, "data", "genomes", "mm10")
homermm10TSSfnm <- file.path(homermm10Annotd, "mm10.tss")

ciceroPosProximalDistalConnfnm <- file.path(
  projd, "data", "snATAC", "mba.whole.sa2subclass.pos.pdc.uniq.coacc.bedpe")
ppdcWithCorfnm <- file.path(
  projd, "data", "snATAC", "mba.whole.sa2subclass.pearson.pos.pdc.alignv1.tsv"
)
eRNAd <- file.path(projd, "13.cicero", "out", "eRNA")
eRNAlogRPKMfnm <- file.path(eRNAd, "eRNA.logRPKM.filteredDAREhn.pbysc.rds")
eRNAfdpath <- file.path(eRNAd, "mat")

mm10chromSizefnm <- file.path(projd, "meta", "mm10.chrom.sizes.lite")

AllenClMetafnm <- file.path(projd, "meta",
  "AIT21_annotation_freeze_081523.tsv")
AllenTFOrdfnm <- file.path(projd, "meta",
  "41586_2023_6812_MOESM9_ESM.txt")

promoterfnm <- file.path("/projects/ps-renlab2/szu/genome",
  "mouseGenePromoter.GRCm38-M25.promoter.1kbAwayTSS.tsv")

# * chromHMM
stateNameOrd <- c(
  "Chr-A", "Chr-B", "Chr-R", "Chr-P", "Chr-O", "Hc-P", "Hc-H", "ND")
chromStateGroupfnm <- file.path(projd, "meta",
  "chromHMM18statesAnnotation.csv")
CREAnnotd <- file.path(projd, "06.ChromHMM", "out", "CREAnnot250429")

#' @export
getChromStateColor <- function() {
  chromStateGroup <- data.table::fread(
    file = file.path(chromStateGroupfnm),
    header = T, sep = ",", data.table = F
  ) |> x => `rownames<-`(x, paste0("raw", x$chromHMM_state))
  t <- chromStateGroup[, c("Group", "Name", "latest color-ZW")] |>
    x => unique(x) |>
    setNames(object = _, nm = c("Group", "Name", "color"))
  stateName2Color <- t$color |> setNames(object = _, nm = t$Name)
  return(stateName2Color)
}

#' @export
loadAllenClMeta <- function() {
  data.table::fread(
    file = AllenClMetafnm, sep = "\t", header = TRUE, data.table = FALSE
  ) |>
    x => `rownames<-`(x, paste0("cl-", x[, 1]))
}

#' Load 260 TFs in order used by Allen for their paper Fig.5.d
#' The raw is directly downloaded from the paper supplemantary file.
#' Then manually remove header and convert to txt file.
#' 
#' @return array of TF gene symbols
#' @export
loadAllenTFOrd <- function(){
  data.table::fread(
    file = AllenTFOrdfnm, header = F, data.table = F)$V1
}

#' @export
getATACSubclassName <- function(rawnm) {
  gsub("^\\d+_", "", x = rawnm) |>
    gsub("  ", " ", x = _) |>
    gsub(" ", "_", x = _) |>
    gsub("/", "-", x = _)
}

#' @export
getPairedTagSubclassName <- function(rawnm) {
  gsub("  ", " ", x = rawnm) |>
    gsub(" ", "_", x = _) |>
    gsub("/", "_", x = _) |>
    gsub("-", "_", x = _)
}

getSMARTSeqSubclassName <- function(rawAllenSubclassname) {}

#' @export
setCondaEnv <- function(envnm = "sa2stable",
                        conda = "/home/szu/mambaforge/bin/conda") {
  reticulate::use_condaenv(
    condaenv = envnm,
    conda = conda,
    required = TRUE
  )
}

#' @export
getmm10ChromSize <- function() {
  data.table::fread(
    file = mm10chromSizefnm,
    header = F, data.table = F
  ) |>
    setNames(object = _, nm = c("chr", "size")) |>
    x => `rownames<-`(x, x$chr)
}

#' @export
loadptscMeta <- function(f = ptscMetafnm) {
  data.table::fread(file = f, sep = ",", data.table = F, header = T) |>
    x => `rownames<-`(x, x$PairedTagName)
}

#' Get cell number record from a typical modality string like
#' "H3K4me1-male125:92:33_female105:79:26"
#'
#' @param x modality string
#' @return integer of the total number of cells
#' @export
getNCellOfModality <- function(x) {
  x1 <- unlist(strsplit(x, "-"))[2]
  xsex <- unlist(strsplit(x1, "_"))
  nmale <- unlist(strsplit(xsex[1], ":"))[1] |>
    x => gsub("male", "", x) |>
    as.integer()
  nfemale <- unlist(strsplit(xsex[2], ":"))[1] |>
    x => gsub("female", "", x) |>
    as.integer()
  return(nmale + nfemale)
}

#' Get PairedTag subclasses' number of cells for different
#' moadlities.
#'
#' Use function of \link{getNCellOfModality} to extract the
#' number information.
#' 
#' @param ptscMeta  see details in \link{loadptscMeta}
#' @return data.frame with columns of
#' PairedTagSubclassName, different modalities and overall DNA numbers
#' @export
getPairedTagSubclass2Nmod <- function(ptscMeta) {
  lapply(ptscMeta$PairedTagName, \(sc) {
    data.frame(
      sc = sc,
      all = getNCellOfModality(ptscMeta[sc, "all"]),
      H3K4me1 = getNCellOfModality(ptscMeta[sc, "H3K4me1"]),
      H3K27ac = getNCellOfModality(ptscMeta[sc, "H3K27ac"]),
      H3K27me3 = getNCellOfModality(ptscMeta[sc, "H3K27me3"]),
      H3K9me3 = getNCellOfModality(ptscMeta[sc, "H3K9me3"])
    )
  }) |> do.call(rbind, args = _) |>
    x => `rownames<-`(x, x$PairedTagName)
}


#' Load default promoter as GenomicRanges object.
#' @return GenomicRanges object
#' @export
loadPromoter <- function(fnm = promoterfnm) {
  r <- data.table::fread(file = fnm, sep = "\t",
    header = F, data.table = F) |>
    x => setNames(
      object = x,
      nm = c("chr", "start", "end", "gene", "strand", "gene_type", "RefSeqRNAid")
    )
  # function from genomicRanges.R
  loadGRfromBed(beds = r, colnms = colnames(r)) 
}
