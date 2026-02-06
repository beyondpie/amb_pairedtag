library(Seurat)
library(tidyverse)
Sys.setenv("_R_USE_PIPEBIND_" = "true")

# * meta
mrs <- c("AMY", "CPU", "HYP", "HIP", "ERC", "NAC", "VTA", "PFC")
projd <- here::here()
workd <- file.path(projd, "03.integration")
tfd <- file.path(workd, "out", "tfneu_vf_region_cca_k5")
outd <- file.path(tfd, "un_coembed")
dir.create(outd, showWarnings = FALSE)
# * functions
load_sumtf_seu <- function(r) {
  readRDS(file.path(tfd, str_glue("sumtfneu_vf_{r}_cca_k5.seu.rds")))
}
load_query_seu <- function(r) {
  file.path(tfd, paste0("tf_", r),
            "query.with.tf-cca-kac5_on-cl.rds") |>
    readRDS(file = _)
}

# * ggplot
mytheme <- theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  plot.title = element_text(colour = "black", hjust = 0.5,
                            size = 15),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_text(colour = "black", size = 12),
  axis.line = element_line(colour = "black"),
  legend.position = "right",
  legend.title = element_text(colour = "black", size = 13),
  legend.text = element_text(colour = "black", size = 12)
)


# * per region, perform cell detector
mr <- "HYP"
sumtf <- load_sumtf_seu(mr)
querySeu <- load_query_seu(mr)
tfscoremat <- querySeu@assays$prediction.score.id@data
tfscoremax <- querySeu$predicted.id.score
sumtf$tfscore <- 0.0
sumtf$tfscore[!sumtf$isRef] <- tfscoremax[
  with(sumtf@meta.data, barcode[!isRef])]


p_tech <- DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                  split.by = "tech", group.by = "tech")
p_tech
p_tech2 <- DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                   group.by = "tech", alpha = 0.3)
p_tech2
p_tfscore <- FeaturePlot(sumtf, reduction = "UMAP", raster = TRUE,
                         features = "tfscore", cells = with(sumtf@meta.data, barcode[!isRef]),
                         cols = c("white", "red"))
p_tech + p_tfscore


p_cs <- CellSelector(DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                             cells = sumtf$barcode[!sumtf$isRef],
                             cols = c("lightblue"), alpha = 0.4))

# p_cs2 <- CellSelector(DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
#                              cells = sumtf$barcode[!sumtf$isRef],
#                              cols = c("blue"), alpha = 0.4))
# p_cs <- c(p_cs, p_cs2)

# p_cs <- read.table(file =  file.path(outd, str_glue("{mr}.barcode.txt")))$V1

sumtf$tech[sumtf$barcode %in% p_cs] <- "unco"
DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                   group.by = "tech", alpha = 0.3) +
  ggtitle(mr)


write.table(x = p_cs, file = file.path(outd, str_glue("{mr}.barcode.txt")),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# * update 

# ** CPU
mr <- "CPU"
sumtf <- load_sumtf_seu(mr)
querySeu <- load_query_seu(mr)
tfscoremat <- querySeu@assays$prediction.score.id@data
tfscoremax <- querySeu$predicted.id.score
sumtf$tfscore <- 0.0
sumtf$tfscore[!sumtf$isRef] <- tfscoremax[
  with(sumtf@meta.data, barcode[!isRef])]


p_tech <- DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                  split.by = "tech", group.by = "tech")
p_tech
p_tech2 <- DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                   group.by = "tech", alpha = 0.3)
p_tech2
p_tfscore <- FeaturePlot(sumtf, reduction = "UMAP", raster = TRUE,
                         features = "tfscore", cells = with(sumtf@meta.data, barcode[!isRef]),
                         cols = c("white", "red"))
p_tech + p_tfscore


p_cs1 <- CellSelector(DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                             cells = sumtf$barcode[!sumtf$isRef],
                             cols = c("grey"), alpha = 0.4))

# p_cs2 <- CellSelector(DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
#                              cells = sumtf$barcode[!sumtf$isRef],
#                              cols = c("blue"), alpha = 0.4))
# p_cs <- c(p_cs, p_cs2)

# p_cs <- read.table(file =  file.path(outd, str_glue("{mr}.barcode.txt")))$V1

sumtf$tech[sumtf$barcode %in% p_cs1] <- "unco1"
DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
        group.by = "tech", alpha = 0.3) +
  ggtitle(mr)
write.table(x = p_cs1, file = file.path(outd, str_glue("{mr}.barcode.1.txt")),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
p_cs2 <- CellSelector(DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                             cells = sumtf$barcode[!sumtf$isRef],
                             cols = c("blue"), alpha = 0.4))
sumtf$tech[sumtf$barcode %in% p_cs2] <- "unco2"
DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
        group.by = "tech", alpha = 0.3) +
  ggtitle(mr)
write.table(x = p_cs2, file = file.path(outd, str_glue("{mr}.barcode.2.txt")),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

p_cs3 <- read.table(file =  file.path(outd, str_glue("{mr}.barcode.txt")))$V1
sumtf$tech[sumtf$barcode %in% p_cs3] <- "unco3"
DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
        group.by = "tech", alpha = 0.3) +
  ggtitle(mr)

# ** HIP
mr <- "HIP"
sumtf <- load_sumtf_seu(mr)
querySeu <- load_query_seu(mr)
tfscoremat <- querySeu@assays$prediction.score.id@data
tfscoremax <- querySeu$predicted.id.score
sumtf$tfscore <- 0.0
sumtf$tfscore[!sumtf$isRef] <- tfscoremax[
  with(sumtf@meta.data, barcode[!isRef])]


p_tech <- DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                  split.by = "tech", group.by = "tech")
p_tech
p_tech2 <- DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                   group.by = "tech", alpha = 0.3)
p_tech2
p_tfscore <- FeaturePlot(sumtf, reduction = "UMAP", raster = TRUE,
                         features = "tfscore", cells = with(sumtf@meta.data, barcode[!isRef]),
                         cols = c("white", "red"))
p_tech + p_tfscore


p_cs1 <- CellSelector(DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                              cells = sumtf$barcode[!sumtf$isRef],
                              cols = c("grey"), alpha = 0.4))
sumtf$tech[sumtf$barcode %in% p_cs1] <- "unco1"
DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
        group.by = "tech", alpha = 0.3) +
  ggtitle(mr)
write.table(x = p_cs1, file = file.path(outd, str_glue("{mr}.barcode.1.txt")),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

p_cs2 <- CellSelector(DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                              cells = sumtf$barcode[!sumtf$isRef],
                              cols = c("blue"), alpha = 0.4))
sumtf$tech[sumtf$barcode %in% p_cs2] <- "unco2"
DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
        group.by = "tech", alpha = 0.3) +
  ggtitle(mr)
write.table(x = p_cs2, file = file.path(outd, str_glue("{mr}.barcode.2.txt")),
            quote = FALSE, row.names = FALSE, col.names = FALSE)


p_cs3 <- CellSelector(DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                              cells = sumtf$barcode[!sumtf$isRef],
                              cols = c("blue"), alpha = 0.4))
sumtf$tech[sumtf$barcode %in% p_cs3] <- "unco3"
DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
        group.by = "tech", alpha = 0.3) +
  ggtitle(mr)
write.table(x = p_cs3, file = file.path(outd, str_glue("{mr}.barcode.3.txt")),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ** VTA
mr <- "VTA"
sumtf <- load_sumtf_seu(mr)
querySeu <- load_query_seu(mr)
tfscoremat <- querySeu@assays$prediction.score.id@data
tfscoremax <- querySeu$predicted.id.score
sumtf$tfscore <- 0.0
sumtf$tfscore[!sumtf$isRef] <- tfscoremax[
  with(sumtf@meta.data, barcode[!isRef])]


p_tech <- DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                  split.by = "tech", group.by = "tech")
p_tech
p_tech2 <- DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                   group.by = "tech", alpha = 0.3)
p_tech2
p_tfscore <- FeaturePlot(sumtf, reduction = "UMAP", raster = TRUE,
                         features = "tfscore", cells = with(sumtf@meta.data, barcode[!isRef]),
                         cols = c("white", "red"))
p_tech + p_tfscore


p_cs1 <- CellSelector(DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                              cells = sumtf$barcode[!sumtf$isRef],
                              cols = c("grey"), alpha = 0.4))
sumtf$tech[sumtf$barcode %in% p_cs1] <- "unco1"
DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
        group.by = "tech", alpha = 0.3) +
  ggtitle(mr)
write.table(x = p_cs1, file = file.path(outd, str_glue("{mr}.barcode.1.txt")),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

p_cs2 <- CellSelector(DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                              cells = sumtf$barcode[!sumtf$isRef],
                              cols = c("blue"), alpha = 0.4))
sumtf$tech[sumtf$barcode %in% p_cs2] <- "unco2"
DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
        group.by = "tech", alpha = 0.3) +
  ggtitle(mr)
write.table(x = p_cs2, file = file.path(outd, str_glue("{mr}.barcode.2.txt")),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ** PFC
mr <- "PFC"
sumtf <- load_sumtf_seu(mr)
querySeu <- load_query_seu(mr)
tfscoremat <- querySeu@assays$prediction.score.id@data
tfscoremax <- querySeu$predicted.id.score
sumtf$tfscore <- 0.0
sumtf$tfscore[!sumtf$isRef] <- tfscoremax[
  with(sumtf@meta.data, barcode[!isRef])]


p_tech <- DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                  split.by = "tech", group.by = "tech")
p_tech
p_tech2 <- DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                   group.by = "tech", alpha = 0.3)
p_tech2
p_tfscore <- FeaturePlot(sumtf, reduction = "UMAP", raster = TRUE,
                         features = "tfscore", cells = with(sumtf@meta.data, barcode[!isRef]),
                         cols = c("white", "red"))
p_tech + p_tfscore


p_cs1 <- CellSelector(DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                              cells = sumtf$barcode[!sumtf$isRef],
                              cols = c("grey"), alpha = 0.4))
sumtf$tech[sumtf$barcode %in% p_cs1] <- "unco1"
DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
        group.by = "tech", alpha = 0.3) +
  ggtitle(mr)
write.table(x = p_cs1, file = file.path(outd, str_glue("{mr}.barcode.1.txt")),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

p_cs2 <- CellSelector(DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                              cells = sumtf$barcode[!sumtf$isRef],
                              cols = c("blue"), alpha = 0.4))
sumtf$tech[sumtf$barcode %in% p_cs2] <- "unco2"
DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
        group.by = "tech", alpha = 0.3) +
  ggtitle(mr)
write.table(x = p_cs2, file = file.path(outd, str_glue("{mr}.barcode.2.txt")),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
p_cs3 <- CellSelector(DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
                              cells = sumtf$barcode[!sumtf$isRef],
                              cols = c("blue"), alpha = 0.4))

sumtf$tech[sumtf$barcode %in% p_cs3] <- "unco3"
DimPlot(sumtf, reduction = "UMAP", raster = TRUE,
        group.by = "tech", alpha = 0.3) +
  ggtitle(mr)
write.table(x = p_cs3, file = file.path(outd, str_glue("{mr}.barcode.3.txt")),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
