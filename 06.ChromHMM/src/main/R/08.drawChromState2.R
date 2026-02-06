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
  "Chr-A", "Chr-B", "Chr-R", "Chr-P", "Chr-O", "Hc-P", "Hc-H", "ND"
)
chrs <- paste0("chr", 1:19)

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

##  this consider the combination of raw chromState and promoter/distal
## we use 100 + raw chromHMMState for promoters
chromState2Name <- data.frame(
  state = c(1:18, 1:18 + 100),
  name = "Ns"
) |>
  mutate(
    .data = _,
    name = chromStateGroup[
      paste0("raw", c(1:18, 1:18)), "Name"
    ]
  ) |>
  x => `rownames<-`(x, paste0("s", x$state))


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
  "Intergenic",
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
) |>
  mutate(
    .data = _,
    name = chromState2Name[paste0("s", groupId), "name"]
  ) |>
  select(.data = _, -groupId) |>
  group_by(.data = _, name) |>
  summarise(
    .data = _,
    across(where(is.numeric), sum),
    .groups = "drop"
  ) |>
  as.data.frame() |>
  x => `rownames<-`(x, x$name) |>
  x => x[stateNameOrd, ]

scsGroupSum <- lapply(scs, \(sc) {
  data.table::fread(
    file = file.path(stated, str_glue("{sc}.stateOrgSum.csv")),
    sep = ",", header = T, data.table = F
  ) |>
    mutate(
      .data = _,
      name = chromState2Name[paste0("s", groupId), "name"]
    ) |>
    select(.data = _, -groupId) |>
    group_by(.data = _, name) |>
    summarise(
      .data = _,
      across(where(is.numeric), sum),
      .groups = "drop"
    ) |>
    as.data.frame() |>
    x => `rownames<-`(x, x$name) |>
    x => x[stateNameOrd, ]
}) |>
  setNames(object = _, nm = scs)

allGroupSum <- data.table::fread(
  file = file.path(stated, "all.stateOrgSum.csv"),
  header = T, sep = ",", data.table = F
) |>
  mutate(
    .data = _,
    name = chromState2Name[paste0("s", groupId), "name"]
  ) |>
  select(.data = _, -groupId) |>
  group_by(.data = _, name) |>
  summarise(
    .data = _,
    across(where(is.numeric), sum),
    .groups = "drop"
  ) |>
  as.data.frame() |>
  x => `rownames<-`(x, x$name) |>
  x => x[stateNameOrd, ]




# * main
## * plot chromHMM emission matrix
r <- allGroupSum |>
  mutate(
    .data = _,
    ATAC = nATAC / totalCount,
    H3K27ac = nH3K27ac / totalCount,
    H3K27me3 = nH3K27me3 / totalCount,
    H3K4me1 = nH3K4me1 / totalCount,
    H3K9me3 = nH3K9me3 / totalCount
  ) |>
  select(
    .data = _,
    -c(
      "nATAC", "nH3K27ac", "nH3K27me3", "nH3K4me1", "nH3K9me3",
      "totalCount"
    )
  )

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
  "chromHMM.state-18-to-8.transition.heatmap.pdf"
), code = {
  print(p)
}, width = 3, height = 5)

## * plot chrom states coverage
stateCounts <- lapply(scsGroupSum, \(i) {
  i[match(stateNameOrd, i$name), "totalCount"]
}) |>
  do.call(cbind, args = _) |>
  x => `rownames<-`(x, stateNameOrd)

## a1 <- lapply(scsGroupSum, funtion(i) {
##   i[match(stateNameOrd, i$name), "totalCount"]
## })
## a2 <- do.call(cbind, args = a1)
## rownames(a2) <- stateNameOrd
## stateCounts <- a2


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
    labels = c(
      "10 kb", "100 kb", "1 Mb", "10 Mb",
      "100 Mb", "1 Gb", "10 Gb"
    ),
    limits = c(4, 10)
  )

withr::with_pdf(
  new = file.path(figd, "chromHMM.state-18-to-8.genomeCoverage.pdf"),
  code = {
    print(p)
  },
  width = 7, height = 6
)

## * plot genome annotations
r <- matrix(0,
  nrow = nrow(allGenomeAnnotSum),
  ncol = length(ordGenomeElements),
  dimnames = list(
    rownames(allGenomeAnnotSum),
    ordGenomeElements
  )
)

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
  "chromHMM.state-18-to-8.genomeAnnot.heatmap.pdf"
), code = {
  print(p)
}, width = 7, height = 5)

# * plot DNA methylation sigals
mC <- "mCG"
# mC <- "mCH"
## scs <- list.files(file.path(workd, "out", "DNAMeth_binSize200"),
##   full.names = F, include.dirs = F, no.. = T
## ) |>
##   y => y[grepl(str_glue("8state2{mC}.b200.avgByCount.csv"), x = y)] |>
##   unique() |>
##   strsplit(x = _, split = "\\.") |>
##   x => lapply(x, \(i){
##     i[1]
##   }) |>
##   unlist()

r <- data.table::fread(
  file = file.path(
    workd, "out", "chromHMM18to8_DNAMe",
    str_glue("subclass.chromHMM18to8.{mC}.avgByCount.csv")
  ), header = T, sep = ",", data.table = F
) |>
  mutate(.data = _, state = factor(state, levels = rev(stateNameOrd)))

p <- ggplot(data = r, aes(x = state, y = .data[[mC]], fill = state)) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_manual(values = stateName2Color) +
  theme(
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black")
  ) +
  theme_minimal(
    base_size = 15
  )

withr::with_pdf(
  new = file.path(figd, str_glue("chromHMM.state18to8.DNA{mC}.pdf")),
  code = {
    print(p)
  },
  width = 7, height = 6
)

# * distinguish non-neuronal cells and neurons
scNN <- unique(r$subclass) |>
  y => y[grepl("_NN$", x = y)]

rNN <- r[r$subclass %in% scNN, ]
rNeu <- r[!(r$subclass %in% scNN), ]


p <- ggplot(data = rNN, aes(x = state, y = .data[[mC]], fill = state)) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_manual(values = stateName2Color) +
  theme(
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 15, colour = "black", hjust = 0.5)
  ) +
  labs(title = str_glue("{mC} in NN subclasses")) +
  theme_minimal(
    base_size = 15
  )

withr::with_pdf(
  new = file.path(figd, str_glue("chromHMM.state18to8.DNA{mC}.NN.pdf")),
  code = {
    print(p)
  },
  width = 7, height = 6
)


p <- ggplot(data = rNeu, aes(x = state, y = .data[[mC]], fill = state)) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_manual(values = stateName2Color) +
  theme(
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 15, colour = "black", hjust = 0.5)
  ) +
  labs(title = str_glue("{mC} in Neu subclasses")) +
  theme_minimal(
    base_size = 15
  )

withr::with_pdf(
  new = file.path(figd, str_glue("chromHMM.state18to8.DNA{mC}.Neu.pdf")),
  code = {
    print(p)
  },
  width = 7, height = 6
)
