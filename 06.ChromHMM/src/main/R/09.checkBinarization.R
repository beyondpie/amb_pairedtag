library(tidyverse)
library(GenomicRanges)
library(tmpRpkg)
options(scipen = 999)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

# * meta
projd <- here::here()
workd <- file.path(projd, "06.ChromHMM")
chromHMMbinaryd <- file.path(workd, "out", "binarizeBed_200")
chrs <- paste0("chr", seq_len(19))
chromSizes <- data.table::fread(
  file = file.path(projd, "meta", "mm10.chrom.sizes.lite"),
  header = F, sep = "\t", data.table = F
) |>
  x => `rownames<-`(x, x[, 1])
outBinaryBedir <- file.path(workd, "out", "orig_chromHMM_binary_bed")
if (!dir.exists(outBinaryBedir)) {
  dir.create(outBinaryBedir)
}

scs <- list.files(chromHMMbinaryd, pattern = "_binary\\.txt") |>
  y => gsub(pattern = "_chr.*binary\\.txt", "", x = y) |>
  unique() |>
  y => sort(x = y, decreasing = F)

# * functions
getChrPos4ChrBin <- function(chr, chrSize, binSize = 200) {
  nbin <- floor(chrSize / binSize)
  startFrom <- cumsum(c(0, rep(binSize, nbin - 1)))
  endTo <- c(startFrom[-1], binSize * nbin)
  ## zero-based, left closed, right open
  data.frame(
    chr = chr,
    startFrom = startFrom,
    endTo = endTo
  )
}

readBinaryMods <- function(sc, chr) {
  data.table::fread(
    file = file.path(chromHMMbinaryd, str_glue("{sc}_{chr}_binary.txt")),
    sep = "\t",
    skip = 1,
    header = T,
    data.table = F
  )
}

transformBin2Bed <- function(vecOfBin, chrbed, chr) {
  r <- rle(vecOfBin)
  pos <- cumsum(r$lengths)

  endTo <- chrbed[pos, 3]
  startFrom <- chrbed[c(1, (pos + 1)[-length(pos)]), 2]

  # filter out zero positions
  index <- r$values > 0
  endTo <- endTo[index]
  startFrom <- startFrom[index]

  data.frame(
    chr = chr,
    startFrom = startFrom,
    endTo = endTo
  )
}

# * chromHMM binary states to bed files
h <- "H3K27me3"
binSize <- 200
chrPos <- lapply(chrs, \(chr) {
  getChrPos4ChrBin(chr, chromSizes[chr, 2], binSize)
}) |> setNames(object = _, nm = chrs)

beds <- lapply(scs, \(sc) {
  lapply(chrs, \(chr) {
    d <- readBinaryMods(sc, chr)
    transformBin2Bed(d[, h], chrPos[[chr]], chr)
  }) |> do.call(what = rbind, args = _)
}) |> setNames(object = _, nm = scs)

## output beds
for (sc in scs) {
  data.table::fwrite(
    x = beds[[sc]],
    file = file.path(outBinaryBedir, str_glue("{sc}_{h}_chromHMM_binary.bed")),
    sep = "\t", col.names = F, row.names = F
  )
}

# * get H3K27me3 peaks to bed with SPM scores
reproPeakd <- file.path(
  projd, "data", "pairedtag_peak",
  "subclass_peak"
)
outd <- file.path(
  projd, "data", "pairedtag_peak",
  "subclass_peakwithSPM_bed"
)
for (sc in scs) {
  fnm <- file.path(reproPeakd, str_glue("{sc}-H3K27me3.bedtoolmerge.peak"))
  if (file.exists(fnm)) {
    r <- data.table::fread(
      file = fnm, sep = "\t", header = T,
      data.table = F
    )
    r2 <- r[, c("chrom", "startFrom", "endTo", "name", "ScorePerMillion")]
    r2$strand <- "+"
    data.table::fwrite(
      x = r2,
      file = file.path(outd, str_glue("{sc}-H3K27me3.peak.SPM.bed")),
      col.names = F, row.names = F, sep = "\t"
    )
  }
}

# * plot hist
q <- 0.9
plotHist <- function(sc, q = 0.9, nbreaks = 100, h = "H3K27me3") {
  r <- data.table::fread(
    file = file.path(reproPeakd, str_glue("{sc}-{h}.bedtoolmerge.peak")),
    header = T, sep = "\t", data.table = F
  )
  r2 <- r$ScorePerMillion
  print(max(r2))
  print(quantile(r2))
  print(quantile(r2, q))
  r2[r2 >= quantile(r2, q)] <- quantile(r2, q)
  r$ScorePerMillion <- r2
  # hist(r2, breaks = nbreaks, main = sc)
  ## use gglot hist
  p <- ggplot(r, aes(x = ScorePerMillion)) +
    geom_histogram(bins = nbreaks, fill = "blue", color = "black") +
    ggtitle(str_glue("{sc} - {h}")) +
    theme_minimal()
  return(p)
}

r1 <- plotHist("319_Astro_TE_NN", q = 0.9, nbreaks = 50)
r4 <- plotHist("319_Astro_TE_NN", q = 0.9, nbreaks = 50, h = "H3K4me1")
r5 <- plotHist("319_Astro_TE_NN", q = 0.9, nbreaks = 50, h = "H3K9me3")


r2 <- plotHist("005_L5_IT_CTX_Glut", nbreaks = 50)



plotHist2 <- function(sc, q = 0.9, nbreaks = 100) {
  r <- data.table::fread(
    file = file.path(reproPeakd, str_glue("{sc}-H3K27me3.bedtoolmerge.peak")),
    header = T, sep = "\t", data.table = F
  )
  r2 <- r$ScorePerMillion
  max(r2)
  quantile(r2)
  quantile(r2, q)

  r <- r[-which(r2 >= quantile(r2, q)), ]
  # hist(r2, breaks = nbreaks, main = sc)
  ## use gglot hist
  p <- ggplot(r, aes(x = ScorePerMillion)) +
    geom_histogram(bins = nbreaks, fill = "blue", color = "black") +
    ggtitle(sc) +
    theme_minimal()
  message(quantile(r[, "ScorePerMillion"], 0.5))
  return(list(fig = p, data = r))
}
r11 <- plotHist2("319_Astro_TE_NN", nbreaks = 50)
r11$fig

r22 <- plotHist2("005_L5_IT_CTX_Glut", nbreaks = 50)
r22$fig

plotHist3 <- function(sc, q = 0.9, nbreaks = 100) {
  r <- data.table::fread(
    file = file.path(reproPeakd, str_glue("{sc}-H3K27ac.BestSPM.peak")),
    header = T, sep = "\t", data.table = F
  )
  r2 <- r$ScorePerMillion
  print(max(r2))
  print(quantile(r2))
  print(quantile(r2, q))
  r2[r2 >= quantile(r2, q)] <- quantile(r2, q)
  r$ScorePerMillion <- r2
  # hist(r2, breaks = nbreaks, main = sc)
  ## use gglot hist
  p <- ggplot(r, aes(x = ScorePerMillion)) +
    geom_histogram(bins = nbreaks, fill = "blue", color = "black") +
    ggtitle(sc) +
    theme_minimal()
  return(p)
}

r3 <- plotHist3("319_Astro_TE_NN", q = 0.9, nbreaks = 50)

## ** filter Histone peaks by SPM.
reproPeakd <- file.path(
  projd, "data", "pairedtag_peak",
  "subclass_peak"
)

outpeakd <- file.path(
  projd, "data", "chromHMM", "subclass_peak_SPMq25"
)

scs <- data.table::fread(
  file = file.path(outpeakd, "cellmarkfiletable.tsv"), sep = "\t",
  header = F, data.table = F
) |>
  x => unique(x[, 1])
chrs <- paste0("chr", seq_len(19))

hs <- c("H3K27ac", "H3K27me3", "H3K4me1", "H3K9me3")
hsuffix <- c("BestSPM", rep("bedtoolmerge", 3)) |>
  setNames(object = _, nm = hs)

for (sc in scs) {
  for (h in names(hsuffix)) {
    fnm <- file.path(reproPeakd, str_glue("{sc}-{h}.{hsuffix[h]}.peak"))
    if (file.exists(fnm)) {
      raw <- data.table::fread(
        file = fnm, header = T,
        sep = "\t", data.table = F
      )
      # keep only chrs peaks
      raw <- raw[raw$chr %in% chrs, ]
      nraw <- nrow(raw)
      s <- raw$ScorePerMillion
      cutoff <- quantile(s, 0.25)
      r <- raw[s >= cutoff, c("chrom", "startFrom", "endTo", "name")]
      nleft <- nrow(r)
      message(paste(
        sc, h, nraw, nleft,
        round(nleft * 100 / nraw, 2), round(cutoff, 2)
      ))
      data.table::fwrite(
        x = r,
        file = file.path(outpeakd, str_glue("{sc}.{h}.peak.bed")),
        sep = "\t", col.names = F, row.names = F
      )
    }
  }
}

# * compare genome coverage for two 18-state models
raw_18model_path <- file.path(workd, "out", "model_bypeak_b200_s18")
raw_emptyState <- "E1"
SPMfilter_18model_path <- file.path(workd, "out", "model_bypeakSPMq25_b200_s18")
SPMfilter_emptyState <- "E18"

SPMfilterState2rawState <- data.frame(
  SPMfilterState = paste0("E", c(18, 16, 3:6, 1, 2, 17, 7:15)),
  rawState = paste0("E", 1:18)
) |> x => `rownames<-`(x, x$SPMfilterState)

scs <- file.path(
  projd, "data",
  "chromHMM", "subclass_peak", "cellmarkfiletable.tsv"
) |>
  data.table::fread(file = _, header = F, sep = "\t", data.table = F) |>
  x => unique(x[, 1]) |>
  sort()

getGenomeCoverage <- function(fromd, sc, state = 18) {
  data.table::fread(
    file = file.path(fromd, str_glue("{sc}_{state}_segments.bed")),
    sep = "\t", header = F, data.table = F
  ) |>
    setNames(object = _, nm = c("chr", "startFrom", "endTo", "state")) |>
    group_by(.data = _, state) |>
    summarise(.data = _, n = sum(endTo - startFrom))
}

rawCovs <- lapply(scs, \(sc) {
  getGenomeCoverage(raw_18model_path, sc)
}) |> setNames(object = _, nm = scs)


SPMfilterCovs <- lapply(scs, \(sc) {
  getGenomeCoverage(SPMfilter_18model_path, sc)
}) |> setNames(object = _, nm = scs)

getTriMat <- function(listOfStats) {
  lapply(names(listOfStats), \(sc) {
    r <- listOfStats[[sc]]
    r$sc <- sc
    return(r)
  }) |> do.call(rbind, args = _)
}

rawTriMat <- getTriMat(rawCovs) |>
  mutate(log10bp = log10(n + 1)) |>
  as.data.frame() |>
  x => `rownames<-`(x, paste0(x$state, x$sc, sep = "-"))

SPMfilterTriMat <- getTriMat(SPMfilterCovs) |>
  mutate(rawState = SPMfilterState2rawState[state, "rawState"]) |>
  select(-state) |>
  mutate(state = rawState) |>
  select(-rawState) |>
  mutate(log10bp = log10(n + 1)) |>
  as.data.frame() |>
  x => `rownames<-`(x, paste0(x$state, x$sc, sep = "-"))

deltaTriMat <- rawTriMat
deltaTriMat$n <- SPMfilterTriMat[
  match(rownames(rawTriMat), rownames(SPMfilterTriMat)), "n"] - rawTriMat$n

deltaTriMat$log10bp <- with(deltaTriMat,
(n / abs(n + 1e-10)) * log10(abs(n+1e-10)))

rawbp <- ggplot(data = rawTriMat, aes(x = state, y = log10bp, fill = state)) +
  geom_boxplot() +
  coord_flip() +
  scale_y_continuous(
    breaks = c(4, 5, 6, 7, 8, 9, 10),
    labels = c(
      "10 kb", "100 kb", "1 Mb", "10 Mb",
      "100 Mb", "1 Gb", "10 Gb"
    ),
    limits = c(4, 10)
  ) +
  xlab("Raw ChromHMM model")

SPMfilterbp <- ggplot(
  data = SPMfilterTriMat,
  aes(x = state, y = log10bp, fill = state)
) +
  geom_boxplot() +
  coord_flip() +
  scale_y_continuous(
    breaks = c(4, 5, 6, 7, 8, 9, 10),
    labels = c(
      "10 kb", "100 kb", "1 Mb", "10 Mb",
      "100 Mb", "1 Gb", "10 Gb"
    ),
    limits = c(4, 10)
  ) +
  xlab("SPMQ25 ChromHMM Model")

deltabp <- ggplot(
  data = deltaTriMat,
  aes(x = state, y = log10bp, fill = state)
) +
  geom_boxplot() +
  coord_flip() +
  xlab("SPMQ25 - Raw")

figd <- file.path(workd, "figure")
withr::with_pdf(
  new = file.path(figd, "chromHMM18-raw-vs-SPMfilter.genomeCoverage.pdf"),
  code = {
    print(ggpubr::ggarrange(rawbp, SPMfilterbp, deltabp, ncol = 3,
      common.legend = T, legend = "none"))
  },
  width = 18, height = 6
)

  

rawCovs <- lapply(scs, \(sc) {
  r <- getGenomeCoverage(raw_18model_path, sc)
  sum(r$n[!(r$state %in% c(raw_emptyState))])
})

SPMfilterCovs <- lapply(scs, \(sc) {
  r <- getGenomeCoverage(SPMfilter_18model_path, sc)
  sum(r$n[!(r$state %in% c(SPMfilter_emptyState))])
})

(avg_rawCovs <- median(unlist(rawCovs)))

(avg_SPMfilterCovs <- median(unlist(SPMfilterCovs)))
