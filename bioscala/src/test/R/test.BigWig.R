library(tidyverse)
library(reticulate)

# * set python env
use_condaenv(condaenv = "sa2",
  conda = "/home/szu/mambaforge/bin/conda")
pybw <- import("pyBigWig")


projd <- here::here()
bwd <- file.path(projd, "data", "ptDNAbam", "bigwig")
sc <- "001_CLA_EPd_CTX_Car3_Glut"
h <- "H3K4me1"
bgf <- file.path(bwd, str_glue("{sc}.{h}.e100.bs100.sm300.bw"))
chr <- "chr2"
startFrom <- 4782300L
endTo <- 4782400L


bgr <- pybw$open(bgf)

bgr$stats(chr, 4780800L, 4782400L)
bgr$stats(chr, 4780800L, 4782400L, "mean")





