library(tidyverse)

# * meta
projd <- here::here()
workd <- file.path(projd, "06.ChromHMM")
conservfnm <- file.path(
  workd, "out",
  "chromHMMStateConservation.csv"
)
outd <- file.path(workd, "out", "varOfState")

# * main
x <- data.table::fread(
  file = conservfnm, sep = ",", header = T, data.table = F
)
scs <- unique(x$sc1)

getVarOfStatMat <- function(bg) {
  y <- x[x$sc1 == bg, ]
  listOfVec <- lapply(seq_len(nrow(y)), \(i) {
    s <- y[i, 3]
    ss <- strsplit(s, split = ";") |> unlist()

    states <- lapply(ss, \(j) {
      strsplit(j, split = ":") |>
        unlist() |>
        x => x[1]
    }) |>
      unlist()
    conservs <- lapply(ss, \(j){
      strsplit(j, split = ":") |>
        unlist() |>
        x => (1.0 - as.numeric(x[2]))
    }) |>
      unlist() |>
      setNames(object = _, nm = states)
  })
  unifnms <- names(listOfVec[[1]])
  m <- lapply(listOfVec, \(i) i[unifnms]) |>
    do.call(rbind, args = _)
  colnames(m) <- unifnms
  rownames(m) <- y$sc2
  return(m)
}

allMat <- lapply(scs, getVarOfStatMat) |>
  setNames(object = _, nm = scs)

saveRDS(allMat, file = file.path(outd, "listOfvarOfState.rds"))

# ** save data from TSCC
# keep using the scs defined above
d <- file.path(workd, "out", "varOfState")

getVarOfStatMat2 <- function(sc) {
  y <- data.table::fread(
    file = file.path(d, str_glue("{sc}.chromHMMStateConservation.Jaccard.csv")),
    sep = ",", header = T, data.table = F
  )
  listOfVec <- lapply(seq_len(nrow(y)), \(i) {
    s <- y[i, 3]
    ss <- strsplit(s, split = ";") |> unlist()

    states <- lapply(ss, \(j) {
      strsplit(j, split = ":") |>
        unlist() |>
        x => x[1]
    }) |>
      unlist()
    conservs <- lapply(ss, \(j){
      strsplit(j, split = ":") |>
        unlist() |>
        x => (1.0 - as.numeric(x[2]))
    }) |>
      unlist() |>
      setNames(object = _, nm = states)
  })
  unifnms <- names(listOfVec[[1]])
  m <- lapply(listOfVec, \(i) i[unifnms]) |>
    do.call(rbind, args = _)
  colnames(m) <- unifnms
  rownames(m) <- y$sc2
  return(m)
}

allMat <- lapply(scs, getVarOfStatMat2) |>
  setNames(object = _, nm = scs)
saveRDS(allMat, file = file.path(outd, "listOfvarOfState.Jaccard.rds"))
