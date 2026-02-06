.libPaths("/opt/homebrew/Caskroom/miniforge/base/lib/R/library")
library(tidyverse)

# * main
projd <- here::here()
workd <- file.path(projd, "17.repressiveMarks")
datafromd <- file.path(workd, "out")
outd <- file.path(workd, "out")

CpGmCGK27me3fnm <- file.path(datafromd, "ptscsCpG.mCG.K27me3.csv")

mCG2K27me3 <- data.table::fread(
  file = CpGmCGK27me3fnm, header = T, sep = ",", data.table = F)
scs <- unique(mCG2K27me3$subclass)

getCor <- function(sc, q = 0.05) {
  m <- mCG2K27me3[mCG2K27me3$subclass == sc, c("mCG", "H3K27me3")]
  x <- m$mCG
  y <- m$H3K27me3
  yind <- (y >= quantile(y, q)) & (y > 0.0) & (y <= quantile(y, 0.99))
  xind <- (x >= quantile(x, q)) & (x > 0.0) & (x <= quantile(x, 0.99))
  ind <- yind & xind
  t <- cor.test(x[ind], y[ind], alternative = "two.sided", method = "pearson")
  data.frame(
    sc = sc,
    pc = t$estimate[1],
    neglot10p =-log10(t$p.value[1])
  )
}

plotCor <- function(sc, q = 0.05) {
  m <- mCG2K27me3[mCG2K27me3$subclass == sc, c("mCG", "H3K27me3")]
  x <- m$mCG
  y <- m$H3K27me3
  yind <- (y >= quantile(y, q)) & (y > 0.0) & (y <= quantile(y, 0.99))
  xind <- (x >= quantile(x, q)) & (x > 0.0) & (x <= quantile(x, 0.99))
  ind <- yind & xind
  plot(x[ind], y[ind])
}

partitionCpG <- function(sc, q = 0.05, low=0.2, high = 0.8) {
  m <- mCG2K27me3[mCG2K27me3$subclass == sc, c("mCG", "H3K27me3")]
  x <- m$mCG
  y <- m$H3K27me3
  yind <- (y >= quantile(y, q)) & (y > 0.0) & (y <= quantile(y, 0.99))
  xind <- (x >= quantile(x, q)) & (x > 0.0) & (x <= quantile(x, 0.99))
  ind <- yind & xind
  x <- x[ind]
  y <- y[ind]
  K27me3Low <- data.frame(
    H3K27me3 = y[x <= 0.2],
    mCGCat = "Low"
  )
  K27me3Mid <- data.frame(
    H3K27me3 = y[ x <= 0.8 & x > 0.2],
    mCGCat = "Medium"
  )
  K27me3High <- data.frame(
    H3K27me3 = y[x > 0.8],
    mCGCat = "High"
  )

  m <- do.call(rbind, list(K27me3Low, K27me3Mid, K27me3High))
  m$mCGCat <- factor(m$mCGCat, levels = c("Low", "Medium", "High"))
  return(m)
}



r <- lapply(scs, getCor) |>
  x => do.call(what = rbind, args = x)

boxplot(r$pc)
boxplot(r$neglog10, ylim = c(2, 100))

sc <- "318_Astro_NT_NN"
plotCor(sc, q = 0.05)

sc <- "327_Oligo_NN"
plotCor(sc, q = 0.05)

sc <- "037_DG_Glut"
plotCor(sc, q = 0.05)

sc <- "053_Sst_Gaba"
plotCor(sc, q = 0.05)
m <- partitionCpG(sc)

(

  p <- ggplot(data = m, aes(x = mCGCat, y = H3K27me3, fill = mCGCat)) +
  geom_boxplot()

)
