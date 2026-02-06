#' 10X gene mtx to Seurat object
#' @export
read10X.seurat <- function(dir_of_10X, minCell = 3, minFeature = 100,
                           project = "10X",
                           mitoPattern = "^mt-") {
  Seurat::Read10X(data.dir = dir_of_10X) |>
    Seurat::CreateSeuratObject(counts = _, project = project,
      min.cells = minCell, min.features = minFeature) |>
    Seurat::PercentageFeatureSet(object = _, pattern = mitoPattern)
}
