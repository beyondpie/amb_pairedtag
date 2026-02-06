library(tidyverse)

projd <- here::here()
rawModeld <- file.path(projd, "data/chromHMM", "model_bypeak_b200_s18")

sc <- "001_CLA_EPd_CTX_Car3_Glut_18"

x <- data.table::fread(
  file = file.path(rawModeld, str_glue("{sc}_dense.bed")),
  header = F,
  sep = "\t",
  skip = 0,
  data.table = F
)


## dense bed model index start from 1 to 18.
# max(x$V4)
# [1] 18

# min(x$V4)
# [1] 1
