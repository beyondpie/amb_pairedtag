library(tidyverse)

projd <- here::here()
workd <- file.path(projd, "17.repressiveMarks")
outd <- file.path(workd, "out", "negCorRepGene")

genes <- data.table::fread(
  file = file.path(outd, "genes.txt"),
  header = F, data.table = F
)[ , 1] |> unlist()

ptscs <- data.table::fread(
  file = file.path(outd, "ptscs.txt"),
  header = F, data.table = F
)[,1] |> unlist()

gExpMat <- data.table::fread(
  file = file.path(outd, "gene2ptscs_RNA.csv"),
  sep = ",",
  header = F, data.table = F
) |>
  x => setNames(object = x, nm = ptscs)


k9Mat <- data.table::fread(
  file = file.path(outd, "gene2ptscs_H3K9me3.csv"),
  sep = ",",
  header = F, data.table = F
) |>
  x => setNames(object = x, nm = ptscs)

k27Mat <- data.table::fread(
  file = file.path(outd, "gene2ptscs_H3K27me3.csv"),
  sep = ",",
  header = F, data.table = F
) |>
  x => setNames(object = x, nm = ptscs)


getCor <- function(gindex, reprsv = "H3K9me3", q = 0.01) {
  x <- gExpMat[gindex, ] |> unlist() |> as.array()
  y <- if(reprsv == "H3K9me3") {
    k9Mat[gindex, ]
  } else {
    k27Mat[gindex, ]
  }
  y <- y |> unlist() |> as.array()
  yind <- (y >= quantile(y, q)) & (y <= quantile(y, 0.99))
  xind <- (x >= quantile(x, q)) & (x <= quantile(x, 0.99))
  ind <- yind & xind
  t <- cor.test(x[ind], y[ind], alternative = "two.sided", method = "pearson")
  data.frame(
    gene = genes[gindex],
    pc = t$estimate[1],
    p = t$p.value[1],
    neglot10p =-log10(t$p.value[1])
  )  
}

corK9onGene <- lapply(seq_along(genes), \(i) {
  getCor(i, reprsv = "H3K9me3", q = 0.01)
}) |>
  do.call(rbind, args = _) |>
  x => x[!is.na(x$pc), ] |>
  x => mutate(.data = x,qval = p.adjust(x$p, method = "BH"))


corK27onGene <- lapply(seq_along(genes), \(i) {
  getCor(i, reprsv = "H3K27me3", q = 0.01)
}) |>
  do.call(rbind, args = _) |>
  x => x[!is.na(x$pc), ] |>
  x => mutate(.data = x,qval = p.adjust(x$p, method = "BH"))

# 2270 genes
negCorK9onGene <- corK9onGene |>
  x => x[!duplicated(x$gene), ] |>
  x => x[ (x$qval <= 0.01) & (x$pc < (-0.1)), ]

# 5829 genes
negCorK27onGene <- corK27onGene |>
  x => x[!duplicated(x$gene), ] |>
  x => x[ (x$qval <= 0.01) & (x$pc < (-0.1)), ]


saveRDS(
  object = list(
    corK9onGene = corK9onGene,
    corK27onGene = corK27onGene,
    negCorK27onGene = negCorK27onGene,
    negCorK9onGene = negCorK9onGene
  ),
  file = file.path(outd, "negCorReprsvOnGenes.rds")
)
