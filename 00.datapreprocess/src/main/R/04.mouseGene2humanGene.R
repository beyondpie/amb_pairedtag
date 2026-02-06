library(tmpRpkg)
proj <- here::here()

Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

# Ref:https://bioinformatics.stackexchange.com/questions/17486/converting-mouse-genes-to-human-genes


homologs <- read.csv(
  "http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",
  sep = "\t")
mouseGeneAnnot <- subset(homologs,
  Common.Organism.Name == "mouse, laboratory") |>
  x => x[!duplicated(x$Symbol), ] |>
  x => `rownames<-`(x, x$Symbol)

humanGeneAnnot <- subset(homologs, Common.Organism.Name == "human") |>
  x => x[!duplicated(x$Symbol), ] |>
  x => `rownames<-`(x, x$Symbol)

ms2hs <- function(gene_list) {
  lapply(gene_list, \(g) {
    c1 <- g %in% rownames(mouseGeneAnnot)
    c2 <- ifelse(c1, any(
      humanGeneAnnot[["DB.Class.Key"]] == mouseGeneAnnot[g, "DB.Class.Key"]), F)
    if (c1 && c2) {
      classKey <- mouseGeneAnnot[g, "DB.Class.Key"]
      humanGenes <- with(humanGeneAnnot, Symbol[DB.Class.Key == classKey])
      hg <- if (length(humanGenes) > 1) {
        humanGenes[1]
      } else {
        humanGenes
      }
      data.frame(
        ms = g,
        hm = hg
      )
    } else {
      NULL
    }
  })
}

msGenes <- mouseGeneAnnot$Symbol
ms2hsdf <- ms2hs(msGenes) |>
  filterNULLfromList() |>
  do.call(rbind, args = _) |>
  setNames(object = _, nm = c("mouse", "human"))
rownames(ms2hsdf) <- ms2hsdf$mouse
data.table::fwrite(ms2hsdf,
  file = file.path(projd, "meta", "mouseGene2humanGeneFromMGI.csv"),
  sep = ",", row.names = F, col.names = T)
