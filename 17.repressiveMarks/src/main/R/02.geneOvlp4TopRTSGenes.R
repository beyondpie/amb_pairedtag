library(tidyverse)
library(tmpRpkg)
# library(inflection)
# library(pspline)
# library(biomaRt)


Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
# * meta
projd <- here::here()
workd <- file.path(projd, "17.repressiveMarks")
rtsd <- file.path(workd, "out")
outd <- file.path(workd, "out")
figd <- file.path(workd, "figure")

# mod <- "H3K27me3"
mod <- "H3K9me3"
qtop <- 0.05

# load Allen's Gene lists
AllenGened <- file.path(projd, "meta", "wmbAllenUpdate")
tfDE <- data.table::fread(file = file.path(AllenGened, "TF534.txt"),
  head = F, data.table = F)$V1
tfMod <- data.table::fread(file = file.path(AllenGened,
  "TF261CoExpModule52.tsv"),
  head = F, sep = "\t", data.table = F) |>
  setNames(object = _, nm = c("mod", "gene", "score")) |>
  x => `rownames<-`(x, x$gene)
gDE <- data.table::fread(file = file.path(AllenGened, "DEG8460.txt"),
  head = F, data.table = F)$V1


# * functions
loadRTS <- function(fnm) {
  data.table::fread(
    file = fnm, sep = "\t", header = F, data.table = F
  ) |> setNames(object = _,
    nm = c("chrom", "startFrom", "endTo", "gene", "rts", "strand"))
}

checkGeneOvlp <- function(gRTS, gTopRTS, gTF, alt = "greater") {
  n11 <- length(intersect(gTopRTS, gTF))
  gNonRTS <- setdiff(gRTS, gTopRTS)
  n21 <- length(intersect(gNonRTS, gTF))
  r <- data.frame(
    "TF" = c(n11, length(gTopRTS) - n11),
    "notTF" = c(n21, length(gNonRTS) - n21),
    row.names = c("RTS", "nonRTS"),
    stringsAsFactors = F
  )
  fisher.test(r, alternative = alt)
}

## * main

# * TODO: qtop selection
# Currently, we use qtop = 0.05 as reported in the literature.
# But this can be chosen based on TF enrichments or other
# criterias.

# * load RTS rankings
rts <- loadRTS(file.path(rtsd,
  str_glue("gene2RTS.{mod}.top{qtop}.all.bed"))) |>
  x => x[!duplicated(x$gene), ] |>
  x => `rownames<-`(x, x$gene)

plot(x = seq_len(nrow(rts)), y = rts$rts)
dy <- abs(diff(rts$rts))
cumdy <- cumsum(dy)
# 727 if 0.9; 1696 if 0.95 for H3K27me3
# 293 if 0.9; 726 if 0.95 for H3K9me3
ng <- min(which(cumdy >= 0.9))
plot(x = seq_len(nrow(rts) - 1), cumdy)
minRTS <- rts$rts[ng]
# now we use minRTS as cutoff here
topRTS <- rts[rts$rts >= minRTS, ] |>
  x => `rownames<-`(x, x$gene)
# save topRTS
data.table::fwrite(topRTS,
file = file.path(outd, str_glue("{mod}.topRTS.gene.bed")),
sep = "\t", col.names = F, row.names = F)

# * check Allen's Gene list Overlapping
checkGeneOvlp(rts$gene, topRTS$gene, tfDE, alt = "greater")
checkGeneOvlp(rts$gene, topRTS$gene, tfMod$gene, alt = "greater")
checkGeneOvlp(rts$gene, topRTS$gene, gDE, alt = "greater")

# * check RTSs reported by the literatures
ms2hs <- data.table::fread(
  file = file.path(projd, "meta", "mouseGene2humanGeneFromMGI.csv"),
  sep = ",", header = T, data.table = F
) |> x => `rownames<-`(x, x$human)

humanRTSd <- file.path(workd, "src/main/resources")
humanRTSEpiMap <- data.table::fread(
  file = file.path(humanRTSd, "human_rts_epimap.txt"),
  header = F, sep = "\t", data.table = F
) |>
  setNames(object = _, nm = c("gene", "rts", "isrts")) |>
  subset(x = _, subset = isrts == "Y") |>
  x => intersect(x$gene, ms2hs$human) |>
  x => ms2hs[x, "mouse"]

humanRTSRoadMap <- data.table::fread(
  file = file.path(humanRTSd, "human_rts_roadmap.txt"),
  header = F, sep = "\t", data.table = F
) |>
  setNames(object = _, nm = c("gene", "rts", "isrts")) |>
  subset(x = _, subset = isrts == "Y") |>
  x => intersect(x$gene, ms2hs$human) |>
  x => ms2hs[x, "mouse"]

checkGeneOvlp(rts$gene, topRTS$gene, humanRTSEpiMap, alt = "greater")
checkGeneOvlp(rts$gene, topRTS$gene, humanRTSRoadMap, alt = "greater")

length(intersect(humanRTSEpiMap, topRTS$gene))
length(intersect(humanRTSRoadMap, topRTS$gene))

# * to add comparision between H3K27me3 vs H3K9me3
