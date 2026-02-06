library(tidyverse)
library(tmpRpkg)

# * meta
projd <- here::here()
ptCREAnnotfnm <- file.path(
  projd, "data", "chromHMM", "allCRE.amb.PairedTag.annot.tsv")
CREHomerAnnotfnm <- file.path(
  projd, "meta", "CREHomerAnnotCleanUp.rds")
ptscMeta <- tmpRpkg::loadptscMeta() |>
  x => x[x$ATAC > 0, ] |>
  x => `rownames<-`(x, x$ATACName)
hTEscfnm <- file.path(projd, "meta", "highTESubclass.csv")

# * main
# ** organize highTE subclasses into a meta file.

## This is to mofidy original single-column file
## after this, the file will be changed to the one
## with two columns.

hTEAllenscs <- data.table::fread(
  file = hTEscfnm,
  header = F, data.table = F) |>
  x => x[, 1]

hTE <- data.frame(
  ATACName = hTEAllenscs,
  PairedTagName = ptscMeta[hTEAllenscs, "PairedTagName"]
)
data.table::fwrite(x = hTE,
  file = hTEscfnm,
  col.names = T, row.names = F, sep = ",")

# ** check the ChrA, ChrO CREs enrichments
