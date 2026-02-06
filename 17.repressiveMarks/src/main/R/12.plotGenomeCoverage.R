library(tidyverse)
library(tmpRpkg)
library(ggpubr)
library(cowplot)

# * meta
projd <- here::here()
workd <- file.path(projd, "17.repressiveMarks")
outd <- file.path(workd, "out")
figd <- file.path(workd, "figure")

# number of cells per subclass
ptscMeta <- tmpRpkg::loadptscMeta()
ptsc2nmod <- tmpRpkg::getPairedTagSubclass2Nmod(ptscMeta)
ptsc2nmod[["Chr-A"]] <- ptsc2nmod$H3K27ac
ptsc2nmod[["Chr-B"]] <- with(ptsc2nmod, H3K4me1 + H3K27me3 + H3K27ac)
ptsc2nmod[["Chr-P"]] <- ptsc2nmod$H3K4me1
ptsc2nmod[["Chr-R"]] <- ptsc2nmod$H3K4me1 + ptsc2nmod$H3K27me3
ptsc2nmod[["Hc-P"]] <- ptsc2nmod$H3K27me3
ptsc2nmod[["Hc-H"]] <- ptsc2nmod$H3K9me3
rownames(ptsc2nmod) <- ptsc2nmod$sc

# * functions
loadGenomeCoverage <- function(tag = "ChromHMM") {
  fnm <- if(tag == "ChromHMM") {
    file.path(outd, "pt.subclassChromHMMGenomeRaioCapture.prcnt.csv")
  } else {
    file.path(outd, "pt.subclassPeakGenomeRaioCapture.prcnt.csv")
  }
  r <- data.table::fread(
    file = fnm, header = T, sep = ",", data.table = F) |>
    x => `rownames<-`(x, x$subclass)
}

plotGenomeCoverage <- function(cvrgMat, ptsc2nmod, useCol = "H3K27me3",
                               needSort = T
                               ) {
  r <- cvrgMat[ , c("subclass", useCol)]

  if(needSort) {
    r1 <- r[!grepl("_NN$", r$subclass), ] |>
      x => x[order(x[[useCol]]), ]
    r2 <- r[grepl("_NN$", r$subclass), ]|>
      x => x[order(x[[useCol]]), ]
    r <- rbind(r1, r2)
  }
  
  d <- ptsc2nmod[r$subclass, c("sc", useCol)]
  r$subclass <- factor(r$subclass, levels = r$subclass)
  d$sc <- factor(d$sc, levels = d$sc)
  d$y <- 1
  d[["log10Count"]] <- log10(d[[useCol]] + 1)
  

  p1 <- ggplot(r, aes(x = .data[["subclass"]], y = .data[[useCol]])) +
    geom_bar(stat = "identity" , fill = "red") +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      panel.border = element_blank()
      )
  p2 <- ggplot(d, aes(x = .data[["sc"]],
    y = .data[["y"]], size = .data[[useCol]])) +
    geom_point() +
    # scale_size_binned(range = c(0.1, 5)) +
    theme(
      axis.text.x = element_text(color = "black", angle = 90, size = 5,
        hjust = 1),
      axis.text.y = element_blank(),
      panel.border = element_blank()

      )
  p <- ggarrange(p1, p2, nrow = 2, ncol = 1, labels = c("A", "B"), align = "h",
    common.legend = T, vjust = 0, heights = c(2, 1), lab)
  return(p)
}


# * main

tag <- "ChromHMM"
cvrgMat <- loadGenomeCoverage(tag = tag)
useCols <- setdiff(colnames(cvrgMat), c("subclass", "Chr-O"))
needSort <- T
for (useCol in useCols) {
  p <- plotGenomeCoverage(cvrgMat, ptsc2nmod, useCol, needSort)
  ggsave(filename = file.path(figd, str_glue("{tag}.{useCol}.pdf")),
    plot = p, width = 15, height = 8)
}

tag <- "Peak"
cvrgMat <- loadGenomeCoverage(tag = tag)
useCols <- setdiff(colnames(cvrgMat), c("subclass", "Chr-O"))
needSort <- T
for (useCol in useCols) {
  p <- plotGenomeCoverage(cvrgMat, ptsc2nmod, useCol, needSort)
  ggsave(filename = file.path(figd, str_glue("{tag}.{useCol}.pdf")),
    plot = p, width = 15, height = 8)
}




