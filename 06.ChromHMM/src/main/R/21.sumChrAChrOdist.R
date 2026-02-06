library(tidyverse)
library(tmpRpkg)

# * meta
projd <- here::here()
workd <- file.path(projd, "06.ChromHMM")
figd <- file.path(workd, "figure")
outd <- file.path(workd, "out")
distd <- file.path(outd, "distOfChrAChrO")

isDAR <- F
suffix <- ifelse(isDAR, "DAR", "all")

# * functions
loadDist <- function(sc, state) {
  fnm <- ifelse(isDAR,
    file.path(distd, str_glue("{sc}.{state}.DAR.dist.csv")),
    file.path(distd, str_glue("{sc}.{state}.dist.csv"))
    )
  data.table::fread(fnm, sep = ",", header = T, data.table = F)
}

# * main
ptscMeta <- tmpRpkg::loadptscMeta()
ptscs <- ptscMeta[ptscMeta$ATAC > 0, "PairedTagName"]

scsChrAdist <- lapply(ptscs, \(sc) {
  loadDist(sc, "Chr-A")
}) |>
  setNames(object = _, nm = ptscs)

ChrAdist <- do.call(rbind, scsChrAdist)

scsChrOdist <- lapply(ptscs, \(sc) {
  loadDist(sc, "Chr-O")
}) |>
  setNames(object = _, nm = ptscs)
ChrOdist <- do.call(rbind, scsChrOdist)

wilcox.test(x = ChrAdist$avgDist, ChrOdist$avgDist, alternative = "less")

r <- rbind(
  data.frame(
    state = "Chr-A",
    dist = ChrAdist$avgDist / 1000
  ),
  data.frame(
    state = "Chr-O",
    dist = ChrOdist$avgDist / 1000
  )
)

(
  p <- ggplot(data = r, aes(x = state, y = dist, fill = state)) +
    geom_violin(color = "black", outlier.shape = NA) +
    geom_boxplot(width = 0.05, color = "black",
      outlier.shape = NA, fill = "white")+
    theme_classic() +
    scale_fill_manual(values = c("#b2182b", "#fee08b")) +
    ggtitle(str_glue("Distances between {suffix} Chr-As or Chr-Os within 1Mbp"))  +
    ylab("Distances (Kbp)") +
    ylim(0, 600)+
    theme(
      axis.text.x = element_text(size = 14, colour = "black"),
      axis.text.y = element_text(size = 14, colour = "black"),
      axis.title.y = element_text(size = 15, colour = "black"),
      plot.title = element_text(size = 15, colour = "black", hjust = 0.5),
      legend.position = "none"
    )
)

saveRDS(r, file = file.path(outd, str_glue("dist.ChrAs.ChrOs.{suffix}.rds")))
ggsave(filename = file.path(figd, str_glue("dist.ChrAs.ChrOs.{suffix}.pdf")),
  plot = p,
  width = 7, height = 7)
