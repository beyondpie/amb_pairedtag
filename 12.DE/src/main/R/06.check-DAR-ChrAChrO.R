# check DAR's ChrA vs ChrO
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
library(tidyverse)
library(tmpRpkg)

# * meta
projd <- here::here()
outd <- file.path(projd, "12.DE", "out")
allCREfnm <- file.path(projd, "data", "chromHMM",
  "allCRE.amb.PairedTag.annot.tsv")

ptscMeta <- tmpRpkg::loadptscMeta() |>
  x => x[x$ATAC > 0, ]
ptscs <- ptscMeta$PairedTagName

# * function
getChrAChrORatio <- function(allCRE, sc) {
  r <- subset(allCRE, subclass == sc)
  rChrA <- subset(r, chromHMMState == "Chr-A")
  rChrO <- subset(r, chromHMMState == "Chr-O")
  data.frame(
    rChrA  = sum(rChrA$isDAR == "DAR") / nrow(rChrA),
    rChrO = sum(rChrO$isDAR == "DAR") / nrow(rChrO),
    rChrADAR = sum(rChrA$isDAR == "DAR") / sum(r$isDAR == "DAR"),
    rChrODAR = sum(rChrO$isDAR == "DAR") / sum(r$isDAR == "DAR"),
    nDChrA = sum(rChrA$isDAR == "DAR"),
    nDChrO = sum(rChrO$isDAR == "DAR"),
    subclass = sc
  )
}

getnAll <- function(s = "all-male159:84:75_female105:96:9") {
  s1 <- strsplit(gsub("all-male|female","", s), "_")[[1]]
  nmale <- as.integer(strsplit(s1[1], ":")[[1]][1])
  nfemale <- as.integer(strsplit(s1[2], ":")[[1]][1])
  nmale + nfemale
}

getnH3K27ac <- function(s = "H3K27ac-male159:84:75_female105:96:9") {
  s1 <- strsplit(gsub("H3K27ac-male|female","", s), "_")[[1]]
  nmale <- as.integer(strsplit(s1[1], ":")[[1]][1])
  nfemale <- as.integer(strsplit(s1[2], ":")[[1]][1])
  nmale + nfemale
  }

# * main
allCRE <- data.table::fread(
  file = allCREfnm, header = T, sep = "\t", data.table = F)

r <- lapply(ptscs, \(sc) {
  getChrAChrORatio(allCRE, sc)}
) |>
  do.call(rbind, args = _)
rownames(r) <- r$subclass

# get H3K27ac cells and total number of cells
# regress the effect of number of cells to ratio
# - plot the ratio vs log of # cells
# - regress out and check the raio box
# - or we can manually filter out some subclasses

nH3K27ac <- vapply(ptscs, \(sc) {
  getnH3K27ac(ptscMeta[sc, "H3K27ac"])
  }, 1L) |>
  x => `names<-`(x, ptscs)
nAll <- vapply(ptscs, \(sc) {
  getnAll(ptscMeta[sc, "all"])
}, 1L) |>
  x => `names<-`(x, ptscs)

d <- data.frame(
  rChrADAR = r$rChrADAR,
  rChrODAR = r$rChrODAR,
  rChrA = r$rChrA,
  rChrO = r$rChrO,
  nH3K27ac = nH3K27ac[r$subclass],
  nAll = nAll[r$subclass],
  nDChrA = r$nDChrA,
  nDChrO = r$nDChrO,
  subclass = r$subclass
)

d2 <- d[d$nH3K27ac >= 1000, ] # nrow(d2) = 73

r2 <- data.frame(
  rChrADAR = d2$rChrADAR,
  rChrODAR = d2$rChrODAR
)
boxplot(r2)

r3 <- data.frame(
  rChrA = d2$rChrA,
  rChrO = d2$rChrO
)
boxplot(r3)

wilcox.test(r2$rChrADAR, r2$rChrODAR, alternative = "greater")
data.table::fwrite(d2, 
  file = file.path(outd, "DAR.ratio.nH3K27ac.csv"),
  sep = ",", row.names = F, col.names = T
)


(
  pChrA_nH3K27ac <- ggplot(d, aes(x = log(nH3K27ac), y = rChrA)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, formula = y ~ x) +
  labs(
    x = "log(nH3K27ac)",
    y = "ChrA ratio",
    title = "ChrA ratio vs log(nH3K27ac)"
  ) +
  theme_bw()
  )

(
  pChrA_nAll <- ggplot(d, aes(x = log(nAll), y = rChrA)) +
    geom_point() +
    geom_smooth(method = "lm", se = T, formula = y ~ x) +
    labs(
      x = "log(nAll)",
      y = "ChrA ratio",
      title = "ChrA ratio vs log(nAll)"
    ) +
    theme_bw()
)

(
  pChrO_nH3K27ac <- ggplot(d, aes(x = log(nH3K27ac), y = rChrO)) +
    geom_point() +
    geom_smooth(method = "lm", se = T, formula = y ~ x) +
    labs(
      x = "log(nH3K27ac)",
      y = "ChrO ratio",
      title = "ChrO ratio vs log(nH3K27ac)"
    ) +
    theme_bw()
  
)

# rChrA vs rChrO
(
  #pChrA_rChrO <- ggplot(d[d$nH3K27ac >= median(d$nH3K27ac), ], aes(x = rChrO, y = rChrA)) +
  pChrA_rChrO <- ggplot(d, aes(x = rChrO, y = rChrA)) +
    geom_point(aes(size = nH3K27ac)) +
    geom_abline(intercept = 0, slope = 1)+
    labs(
      x = "rChrO",
      y = "rChrA",
      title = "rChrA vs rChrO"
    ) +
    theme_bw() +
    scale_size_continuous(
      range = c(1, 10),
      name = "nH3K27ac"
    )
)

# rChrA vs rChrO
(
  pChrA_rChrO <- ggplot(d[d$nH3K27ac >= median(d$nH3K27ac), ], aes(x = nDChrO, y = nDChrA)) +
  # pChrA_rChrO <- ggplot(d, aes(x = nDChrO, y = nDChrA)) +
    geom_point(aes(size = nH3K27ac)) +
    geom_abline(intercept = 0, slope = 1)+
    labs(
      x = "rChrO",
      y = "rChrA",
      title = "nDChrA vs nDChrO"
    ) +
    theme_bw() +
    scale_size_continuous(
      range = c(1, 10),
      name = "nH3K27ac"
    )
)

dd <- d[d$nH3K27ac >= median(d$nH3K27ac), ]
a <- dd$nDChrA / dd$nDChrO
boxplot(log(a))
# fit linear regression
lmChrA <- lm( rChrA ~ log(nH3K27ac), data = d)
summary(lmChrA)

d$mrChrA <- resid(lmChrA) + lmChrA$coefficients[1]

# mrChrA vs rChrO
(
  pmrChrA_rChrO <- ggplot(d, aes(x = rChrO, y = mrChrA)) +
    geom_point(aes(size = nH3K27ac)) +
    geom_abline(intercept = 0, slope = 1)+
    labs(
      x = "rChrO",
      y = "mrChrA",
      title = "mrChrA vs rChrO"
    ) +
    theme_bw() +
    scale_size_continuous(
      range = c(1, 10),
      name = "nH3K27ac"
    )
)

# mrChrA vs rChrA
(
  pmrChrA_rChrA <- ggplot(d, aes(x = rChrA, y = mrChrA)) +
    geom_point(aes(size = nH3K27ac)) +
    geom_abline(intercept = 0, slope = 1)+
    labs(
      x = "rChrA",
      y = "mrChrA",
      title = "mrChrA vs rChrA"
    ) +
    theme_bw() +
    scale_size_continuous(
      range = c(1, 10),
      name = "nH3K27ac"
    )
)

# rChrADAR vs rChrODAR
(
  pChrADAR_rChrODAR <- ggplot(d, aes(x = rChrODAR, y = rChrADAR)) +
    geom_point(aes(size = nH3K27ac)) +
    geom_abline(intercept = 0, slope = 1)+
    labs(
      x = "rChrODAR",
      y = "rChrADAR",
      title = "rChrADAR vs rChrODAR"
    ) +
    theme_bw() +
    scale_size_continuous(
      range = c(1, 10),
      name = "nH3K27ac"
    )
)

# save the plot
ggsave(
  filename = file.path(outd, "rChrA_vs_rChrO_nH3K27ac.pdf"),
  plot = pChrA_rChrO,
  width = 6, height = 4
)
