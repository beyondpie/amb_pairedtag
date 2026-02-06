suppressMessages({
  suppressWarnings({
    library(tidyverse)
    library(GenomicRanges)
    library(ComplexHeatmap)
    library(tmpRpkg)
  })
})

# * meta
projd <- here::here()
workd <- file.path(projd, "06.ChromHMM")
figd <- file.path(workd, "figure")
stated <- file.path(workd, "out", "m-bpeak-s18_pd-obs")
genomeAnnotd <- file.path(workd, "out", "m-bpeak-s18_homer_bannot")
mods <- c("ATAC", "H3K27ac", "H3K4me1", "H3K27me3", "H3K9me3")
stateNameOrd <- c(
  "Pr-A", "Pr-B", "Pr-R", "Pr-Oac", "Pr-OA",
  "En-Sd", "En-Bd", "En-Rd", "En-Oacd",
  "Hc-P", "Hc-PH", "Hc-H", "Ns", "Unk", "Cr-OAd"
)
chromStateGroup <- data.table::fread(
  file = file.path(
    workd, "out", "model_bypeak_b200_s18",
    "org.chromHMM18states.csv"
  ),
  header = T, sep = ",", data.table = F
) |>
  x => x[, c("Group", "Name", "color")] |>
  x => unique(x) |>
  x => `rownames<-`(x, paste0("gid", x$Group))
stateName2Color <- chromStateGroup$color |>
  setNames(object = _, nm = chromStateGroup$Name)

scs <- list.files(
  path = file.path(stated),
  full.names = F, no.. = T
) |>
  y => gsub("_chr\\d+.*", "", y) |>
  y => gsub("\\.stateOrgSum.csv", "", y) |>
  y => y[grepl("\\d+", y)] |>
  x => unique(x)

allGroupSum <- data.table::fread(
  file = file.path(stated, "all.stateOrgSum.csv"),
  header = T, sep = ",", data.table = F
) |> x => `rownames<-`(x, paste0("gid", x$groupId))

ordGenomeElements <- c(
  "CpG",
  "5UTR",
  "TSS",
  "Exon",
  "Intron",
  "TES",
  "3UTR",
  "GAP",
  "LINE",
  "SINE",
  "LTR",
  "SimpleRepeat",
  "Satellite",
  "rRNA",
  "tRNA"
)

allGenomeAnnotSum <- data.table::fread(
  file = file.path(genomeAnnotd, "chromHMM.groupId2HomerGenomeAnnot.csv"),
  sep = ",",
  head = T,
  data.table = F
) |> x => `rownames<-`(x,
  chromStateGroup[paste0("gid", x$groupId), "Name"]) |>
  x => x[stateNameOrd, ]


scsGroupSum <- lapply(scs, \(sc) {
  data.table::fread(
    file = file.path(stated, str_glue("{sc}.stateOrgSum.csv")),
    sep = ",", header = T, data.table = F
  ) |> x => `rownames<-`(x, paste0("gid", x$groupId))
}) |>
  setNames(object = _, nm = scs)



# * main
## * ggplot config
## mytheme <- theme(
##   panel.border = element_blank(),
##   panel.grid.major = element_blank(),
##   panel.grid.minor = element_blank(),
##   panel.background = element_blank(),
##   plot.title = element_text(colour = "black", hjust = 0.5,
##     size = 12),
##   axis.text = element_blank(),
##   axis.ticks = element_blank(),
##   axis.title = element_blank(),
##   axis.line = element_line(colour = "black"),
##   legend.position = "none"
## )

## * panel: chromHMM emission matrix.
r <- allGroupSum |>
  mutate(
    .data = _,
    ATAC = nATAC / totalCount,
    H3K27ac = nH3K27ac / totalCount,
    H3K27me3 = nH3K27me3 / totalCount,
    H3K4me1 = nH3K4me1 / totalCount,
    H3K9me3 = nH3K9me3 / totalCount,
    name = chromStateGroup[paste0("gid", groupId), "Name"]
  ) |>
  select(
    .data = _,
    -c(
      "nATAC", "nH3K27ac", "nH3K27me3", "nH3K4me1", "nH3K9me3",
      "totalCount", "groupId"
    )
  ) |>
  x => `rownames<-`(x, x$name) |>
  x => x[stateNameOrd, ]

probColfn <- circlize::colorRamp2(c(0, 1), c("white", "blue"))
p <- ComplexHeatmap::Heatmap(
  matrix = as.matrix(r[, mods]),
  name = "Prob",
  col = probColfn,
  show_row_names = T,
  show_column_names = T,
  cluster_rows = F,
  cluster_columns = F,
  row_names_gp = gpar(fontsize = 15),
  column_names_gp = gpar(fontsize = 15)
)
withr::with_pdf(new = file.path(
  figd,
  "chromHMM.state-18-to-15.transition.heatmap.pdf"
), code = {
  print(p)
}, width = 3, height = 5)

## * panel: boxplots of chrom states coverage in the genome
stateCounts <- lapply(scsGroupSum, \(i) {
  mutate(
    .data = i,
    name = chromStateGroup[paste0("gid", i$groupId), "Name"]
  ) |>
    x => x[match(stateNameOrd, x$name), "totalCount"]
}) |>
  do.call(cbind, args = _) |>
  x => `rownames<-`(x, stateNameOrd)

longStateCount <- data.frame(
  name = rep(rownames(stateCounts), times = ncol(stateCounts)),
  subclass = rep(colnames(stateCounts), each = nrow(stateCounts)),
  log10bp = log10(as.vector(as.matrix(stateCounts + 1))) + log10(200)
) |>
  mutate(.data = _, name = factor(name, levels = rev(stateNameOrd)))

p <- ggplot(data = longStateCount, aes(x = name, y = log10bp, fill = name)) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_manual(values = stateName2Color) +
  scale_y_continuous(
    breaks = c(4, 5, 6, 7, 8, 9, 10),
    labels = c("10 kb", "100 kb", "1 Mb", "10 Mb",
      "100 Mb", "1 Gb", "10 Gb"),
    limits = c(4, 10)
  )

withr::with_pdf(
  new = file.path(figd, "chromHMM.state-18-to-15.genomeCoverage.pdf"),
  code = {
    print(p)
  },
  width = 7, height = 6
)


## * panel: chromHMM to genomic elements
r <- matrix(0, nrow = nrow(allGenomeAnnotSum),
  ncol = length(ordGenomeElements),
  dimnames = list(
    rownames(allGenomeAnnotSum),
    ordGenomeElements
  ))

ntotal <- allGenomeAnnotSum$nTotal + 0.0001

for (j in ordGenomeElements) {
  r[, j] <- allGenomeAnnotSum[, j] / ntotal
}


probColfn <- circlize::colorRamp2(c(0, 1), c("white", "red"))
p <- ComplexHeatmap::Heatmap(
  matrix = r,
  name = "Prob",
  col = probColfn,
  show_row_names = T,
  show_column_names = T,
  cluster_rows = F,
  cluster_columns = F,
  row_names_gp = gpar(fontsize = 15),
  column_names_gp = gpar(fontsize = 15)
)
withr::with_pdf(new = file.path(
  figd,
  "chromHMM.state-18-to-15.genomeAnnot.heatmap.pdf"
), code = {
  print(p)
}, width = 7, height = 5)
