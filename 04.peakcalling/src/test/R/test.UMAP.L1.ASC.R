# check clustering L1 and L2 for ASC group
library(tidyverse)

projd <- "/tscc/projects/ps-renlab2/szu/projects/amb_pairedtag"
cellMeta <- file.path(projd, "meta", 
                      "pairedtag.cell.meta.all.240626.csv") |>
  data.table::fread(file = _, sep = ",", header = TRUE, data.table = FALSE)

dCellMeta <- cellMeta[sample(x =1: nrow(cellMeta), size = 30000, replace = F), ]
dCellMeta$cluster.l1.id <- factor(dCellMeta$cluster.l1.id)

newL2Id <- levels(dCellMeta$cluster.l1.id)[dCellMeta$cluster.l1.id]
index <- newL2Id != 11
newL2Id[newL2Id != 11] <- -1
index <- newL2Id > 0
newL2Id[index] <- dCellMeta$cluster.l2.id[index]

dCellMeta$cluster.l2.id <- factor(newL2Id)

pL1 <- ggplot(data = dCellMeta, aes(x = cluster.l1.umap.x, y = cluster.l1.umap.y, color = cluster.l1.id)) +
  geom_point()
pL2 <- ggplot(data = dCellMeta, aes(x = cluster.l1.umap.x, y = cluster.l1.umap.y, color = cluster.l2.id)) +
  geom_point(alpha = 0.2, size = 0.1) + 
  scale_color_brewer(palette = "Set1")
