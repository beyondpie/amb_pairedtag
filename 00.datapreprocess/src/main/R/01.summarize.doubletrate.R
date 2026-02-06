library(tidyverse)

projdir <- here::here()
workdir <- file.path(projdir, "00.datapreprocess")

doubletsumfnm <- file.path(workdir,
  "doublet_rate_calculation",
  "mm10_hg38_doublet_20240213.xlsx")

doubletsums <- readxl::excel_sheets(doubletsumfnm) |>
  set_names() |>
  purrr::map(readxl::read_xlsx, path = doubletsumfnm)

doubletRates <- vapply(doubletsums, \(x) {
  as.numeric(x[8, 2])
}, 0.0) |>
  as_tibble(rownames = "exp") |>
  as.data.frame()

dr.exp14 <- data.frame(exp = "Exp14",
  value = with(doubletRates, mean(value[!(exp %in% c("Exp3"))]))
)

dr.est <- rbind(doubletRates, dr.exp14)
write.table(x = dr.est, file = file.path(
  projdir, "meta", "doublet.rate.est.csv"),
  quote = FALSE, sep = ",", row.names = FALSE,
  col.names = TRUE)
