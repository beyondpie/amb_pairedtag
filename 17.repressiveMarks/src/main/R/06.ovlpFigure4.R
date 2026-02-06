library(tidyverse)
library(tmpRpkg)

# * meta
projd <- here::here()

# * main

# * load H3K27me3-marked genes
topK27me3RTS <- data.table::fread(
  file = file.path(
    projd, "17.repressiveMarks",
    "out",
    str_glue("H3K27me3.topRTS.gene.bed")
  ),
  sep = "\t", header = F, data.table = F
) |>
  setNames(object = _,
    nm = c("chrom", "startFrom", "endTo",
"gene", "rts", "strand")) |>
  x => x[!duplicated(x$gene), ] |>
  x => `rownames<-`(x, x$gene)

# * load H3K27ac-marked genes
tfmotif <- readRDS(file.path(projd, "99.figures",
   "out/fig4", "panelE.sc2motif.enrich.rds"))

tfExp <- readRDS(file.path(
  projd,
  "99.figures", "out/fig4", "panelE.sc2motif.geneExp.rds"
))
K27acTFs <- colnames(tfExp)

# * overlap between them
## "Ets1"  "Fli1"  "Pax5"  "Ebf1"  "Esrrb" "Cux2"  "Lhx2"
## "Irf8"  "Lef1"  "Otx2"
## "Sox2"  "Olig2" "Tead1" "Meis1" "Isl1"  "Zic2"  "Lhx6"
## "Pax8"  "Lhx9"  "Creb5"
## [21] "Spi1"
g <- intersect(topK27me3RTS$gene, K27acTFs)

t <- tfExp[!(rownames(tfExp) %in%
               c("026_NLOT_Rho_Glut",
                 "027_L6b_EPd_Glut",
                 "028_L6b_CT_ENT_Glut",
                 "029_L6b_CTX_Glut",
                 "030_L6_CT_CTX_Glut",
                 "031_CT_SUB_Glut",
                 "032_L5_NP_CTX_Glut",
                 "033_NP_SUB_Glut",
                 "034_NP_PPP_Glut"
                 )), colnames(tfExp) %in% g]


low.val.col <- -0.5
high.val.col <- 1.5

legend_labels <- c(round(low.val.col, 1),
  round(high.val.col, 1))

col_fun <- circlize::colorRamp2(
  seq(low.val.col, high.val.col, length = 60),
  viridis::viridis(60)
)

hmTFExp <- ComplexHeatmap::Heatmap(
  matrix = t,
  col = col_fun,
  cluster_columns = F,
  cluster_rows = F,
  show_row_names = T,
  row_names_gp = grid::gpar(fontsize = 9),
  column_names_gp = grid::gpar(fontsize = 10),
  show_column_names = T,
  top_annotation = NULL,
  left_annotation = NULL,
  use_raster = T,
  show_heatmap_legend = T,
  na_col = "#440154FF",
  heatmap_legend_param = list(
    title = "scaled logCPM",
    at = c(low.val.col, high.val.col),
    labels = legend_labels,
    direction = "horizontal"
  )
)
hmTFExp

# * choose representative subclasses
scs <- c(
  "009_L2_3_IT_PIR_ENTl_Glut",
  "016_CA1_ProS_Glut",
  "017_CA3_Glut",
  "025_CA2_FC_IG_Glut"
  "036_HPF_CR_Glut",
  "037_DG_Glut",
  "044_OB_Dopa_Gaba",
  "049_Lamp5_Gaba",
  "061_STR_D1_Gaba",
  "062_STR_D2_Gaba",
  "081_ACB_BST_FS_D1_Gaba",
  "086_MPO_ADP_Lhx8_Gaba",
  "097_PVHd_SBPV_Six3_Prox1_Gaba",
  "124_MPN_MPO_PVpo_Hmx2_Glut",
  "142_HY_Gnrh1_Glut",
  "156_MB_ant_ve_Dmrta2_Glut",
  "180_SCiw_Pitx2_Glut",
  "216_MB_MY_Tph2_Glut_Sero",
  "318_Astro_NT_NN",
  "319_Astro_TE_NN",
  "326_OPC_NN",
  "327_Oligo_NN",
  "330_VLMC_NN",
  "334_Microglia_NN"
)
