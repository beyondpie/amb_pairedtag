## This is an old version AnnotateCREs
## We can now ignore this file

library(tidyverse)
library(viridis)
library(hrbrthemes)
library(tmpRpkg)


# * meta
projd <- here::here()
workd <- file.path(projd, "13.cicero")
annotd <- file.path(workd, "out", "ATACPeakChromHMMAnnot")
scs <- list.files(annotd) |>
  strsplit(x = _, split = "\\.") |>
  vapply(X = _, \(i) i[1], "a") |>
  sort() |>
  unique()
exts <- c(250, 500, 750, 1000)
figd <- file.path(workd, "figure")
stateNameOrd <- c(
  "Chr-A", "Chr-B", "Chr-R", "Chr-P", "Chr-O", "Hc-P", "Hc-H", "ND"
)
chromStateGroup <- data.table::fread(
  file = file.path(
    projd, "meta", "chromHMM18statesAnnotation.csv"
  ),
  header = T, sep = ",", data.table = F
) |> x => `rownames<-`(x, paste0("raw", x$chromHMM_state))

t <- chromStateGroup[, c("Group", "Name", "latest color-ZW")] |>
  x => unique(x) |>
  setNames(object = _, nm = c("Group", "Name", "color"))
stateName2Color <- t$color |> setNames(object = _, nm = t$Name)
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

homerd <- file.path("/projects/ps-renlab2/szu/softwares/homer")
homermm10Annotd <- file.path(homerd, "data/genomes", "mm10")
homermm10TSSfnm <- file.path(homermm10Annotd, "mm10.tss")
homerTSS <- readHomerTSS(homermm10TSSfnm, bufferSize = 2000) |>
  x => x[x$chr %in% chromSizes$chr, ] |>
  x => loadGRfromBed(beds = x, colnms = colnames(x))

# follow the definition of ATAC-seq paper
homerPromoter <- GenomicRanges::promoters(
  x = homerTSS, upstream = 1000, downstream = 500
)

# * functions
readPeakAnnot <- function(sc, ext) {
  r <- data.table::fread(
    file = file.path(
      annotd,
      str_glue("{sc}.getRepState.ATAC.peakSummit.ext{ext}.csv")
    ),
    header = T,
    sep = ",", data.table = F
  )
  a <- strsplit(r[, 2], split = ":") |>
    do.call(rbind, args = _) |>
    as.data.frame() |>
    setNames(object = _, nm = c("slen", "state")) |>
    mutate(.data = _, slen = as.integer(slen))
  pGR <- transformRegion2GR(region = r$CRE) |>
    resizeIgnoreStrand(gr = _, fix = "center", extsize = 250)
  hits <- GenomicRanges::findOverlaps(query = pGR, subject = homerPromoter)
  proxIndex <- S4Vectors::queryHits(hits) |>
    sort() |>
    unique()
  pd <- rep("distal", length(pGR))
  pd[proxIndex] <- "proximal"
  w <- data.frame(
    CRE = r$CRE,
    state = a$state,
    slen = a$slen,
    pd = pd
  )
  return(w)
}

statAnnot <- function(listOfAnnots, usePD = F) {
  scs <- names(listOfAnnots)
  lapply(scs, \(sc) {
    a <- listOfAnnots[[sc]]
    if (usePD) {
      table(with(a, paste(state, pd, sep = ":"))) |>
        as.data.frame(x = _, stringsAsFactors = F) |>
        setNames(object = _, nm = c("statepd", "freqpd")) |>
        mutate(.data = _, sc = sc)
    } else {
      table(a$state) |>
        as.data.frame(x = _, stringsAsFactors = F) |>
        setNames(object = _, nm = c("state", "freq")) |>
        mutate(.data = _, sc = sc)
    }
  }) |>
    do.call(rbind, args = _)
}

# * main
# ** load all the annotations
p_e250 <- lapply(scs, \(sc) {
  readPeakAnnot(sc, ext = 250)
}) |> setNames(object = _, nm = scs)

p_e500 <- lapply(scs, \(sc) {
  readPeakAnnot(sc, ext = 500)
}) |> setNames(object = _, nm = scs)


p_e750 <- lapply(scs, \(sc) {
  readPeakAnnot(sc, ext = 750)
}) |> setNames(object = _, nm = scs)

p_e1000 <- lapply(scs, \(sc) {
  readPeakAnnot(sc, ext = 1000)
}) |> setNames(object = _, nm = scs)

# ** stat different annotations
statAnnot_e250 <- statAnnot(p_e250)
statAnnot_e500 <- statAnnot(p_e500)
statAnnot_e750 <- statAnnot(p_e750)
statAnnot_e1000 <- statAnnot(p_e1000)

statAnnots <- list(
  ext250 = statAnnot_e250,
  ext500 = statAnnot_e500,
  ext750 = statAnnot_e750,
  ext1000 = statAnnot_e1000
)

sumStatAnnot <- lapply(names(statAnnots), \(ext){
  statAnnots[[ext]] |>
    group_by(.data = _, state) |>
    summarise(
      .data = _,
      total_freq = sum(freq)
    ) |>
    setNames(object = _, nm = c("state", "freq")) |>
    as.data.frame() |>
    mutate(
      .data = _,
      ext = ext,
      percnt = freq * 100 / sum(freq)
    )
}) |>
  do.call(rbind, args = _) |>
  mutate(
    .data = _,
    ext = factor(ext, levels = paste0("ext", exts)),
    state = factor(state,
      levels = c(
        "Chr-A", "Chr-B", "Chr-R", "Chr-P", "Chr-O", "Hc-P", "Hc-H",
        "ND"
      )
    )
  )

p <- ggplot(data = sumStatAnnot, aes(x = ext, y = percnt, fill = state)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = stateName2Color) +
  ggtitle("ATAC peaks' annotations using ChromHMM State") +
  theme(
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 15, colour = "black", hjust = 0.5)
  ) +
  theme_minimal(
    base_size = 15
  )

withr::with_pdf(
  new = file.path(figd, str_glue("ATACPeak.sumStateCount4exts.pdf")),
  code = {
    print(p)
  },
  width = 7, height = 6
)

# ** check distal and promoter
statAnnotpd_e250 <- statAnnot(p_e250, T)
statAnnotpd_e500 <- statAnnot(p_e500, T)
statAnnotpd_e750 <- statAnnot(p_e750, T)
statAnnotpd_e1000 <- statAnnot(p_e1000, T)

statAnnotspd <- list(
  ext250 = statAnnotpd_e250,
  ext500 = statAnnotpd_e500,
  ext750 = statAnnotpd_e750,
  ext1000 = statAnnotpd_e1000
)

sumStatAnnotpd <- lapply(names(statAnnotspd), \(ext){
  statAnnotspd[[ext]] |>
    group_by(.data = _, statepd) |>
    summarise(
      .data = _,
      total_freq = sum(freqpd)
    ) |>
    setNames(object = _, nm = c("state", "freq")) |>
    as.data.frame() |>
    mutate(
      .data = _,
      ext = ext,
      percnt = freq * 100 / sum(freq)
    )
}) |>
  do.call(rbind, args = _) |>
  mutate(
    .data = _,
    ext = factor(ext, levels = paste0("ext", exts)),
  )

p <- ggplot(data = sumStatAnnotpd, aes(x = ext, y = percnt, fill = state)) +
  geom_bar(position = "stack", stat = "identity") +
  ggtitle("ATAC peaks' annotations using ChromHMM State") +
  theme(
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 15, colour = "black", hjust = 0.5)
  ) +
  theme_minimal(
    base_size = 15
  )

withr::with_pdf(
  new = file.path(figd, str_glue("ATACPeak.sumStateCount4exts.pd.pdf")),
  code = {
    print(p)
  },
  width = 7, height = 6
)

# ** check total ATAC peaks
atacChrADistal_ext750 <- lapply(p_e750, \(d){
  with(d, CRE[state == "Chr-A" & pd == "distal"])
}) |>
  unlist() |>
  sort() |>
  unique()

atacChrRDistal_ext750 <- lapply(p_e750, \(d){
  with(d, CRE[state == "Chr-R" & pd == "distal"])
}) |>
  unlist() |>
  sort() |>
  unique()

atacChrADistal_ext1000 <- lapply(p_e1000, \(d){
  with(d, CRE[state == "Chr-A" & pd == "distal"])
}) |>
  unlist() |>
  sort() |>
  unique()

atacChrRDistal_ext1000 <- lapply(p_e1000, \(d){
  with(d, CRE[state == "Chr-R" & pd == "distal"])
}) |>
  unlist() |>
  sort() |>
  unique()


# * now plot the annotations from Summit
pSummits <- lapply(scs, \(sc) {
  readSummitAnnot(
    file.path(annotd, str_glue("{sc}.getRepStateOnSummit.ATAC.csv")),
    homerPromoter)
}) |> setNames(object = _, nm = scs)

# statAnnot_pSummit <- list(summit = statAnnot(pSummits, usePD = F))
statAnnot_pSummit <- list(summit = statAnnot(pSummits, usePD = T))

sumStatAnnot <- lapply(names(statAnnot_pSummit), \(ext){
  statAnnot_pSummit[[ext]] |>
    group_by(.data = _, state) |>
    # group_by(.data = _, statepd) |>
    summarise(
      .data = _,
      # total_freq = sum(freqpd)
      total_freq = sum(freq)
    ) |>
    setNames(object = _, nm = c("state", "freq")) |>
    as.data.frame() |>
    mutate(
      .data = _,
      ext = ext,
      percnt = freq * 100 / sum(freq)
    )
}) |>
  do.call(rbind, args = _) |>
  mutate(
    .data = _,
    state = factor(state,
      levels = c(
        "Chr-A", "Chr-B", "Chr-R", "Chr-P", "Chr-O", "Hc-P", "Hc-H",
        "ND"
      )
    )
  )

p <- ggplot(data = sumStatAnnot, aes(x = ext, y = percnt, fill = state)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = stateName2Color) +
  #ggtitle("ATAC peaks' annotations using ChromHMM State") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 15, colour = "black", hjust = 0.5)
  ) +
  theme_minimal(
    base_size = 15
  )

withr::with_pdf(
  new = file.path(figd, str_glue("ATACPeak.sumStateCountBySummit.pdf")),
  # new = file.path(figd,
  #   str_glue("ATACPeak.sumStateCountBySummitProximalDistal.pdf")),
  code = {
    print(p)
  },
  width = 3, height = 6
)

# * prepare Chr-A / Chr-R ATACPeak
# - each subclass
# - combinations of all
# - double check using track data

allPeaks <- data.table::fread(
  file = file.path(projd, "data", "snATAC", "wmb.snATAC.peak.srt.bed"),
  header = F,
  sep = "\t",
  data.table = F
) |> setNames(object = _, nm = c("chr", "start", "end", "name")) |>
  loadGRfromBed(beds = _)


