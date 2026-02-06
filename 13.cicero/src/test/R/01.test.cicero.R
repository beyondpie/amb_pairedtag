suppressMessages({
  suppressWarnings({
    library(cicero)
    library(tidyverse)
    library(reticulate)
    library(Matrix)
    library(tmpRpkg)
    library(future.apply)
  })
})

# * meta
setCondaEnv(envnm = "sa2stable")
future::plan(strategy = multisession, workers = 10)
# use convert as FALSE,
# otherwise hard to visualize the info of ann$X
ad <- reticulate::import(module = "anndata", convert = FALSE)
projd <- here::here()

# * test csr or csc sparse matrix loading to R
scipy <- reticulate::import(module = "scipy", convert = FALSE)
np <- reticulate::import(module = "numpy", convert = FALSE)

id_row <- c(0L, 0L, 1L, 2L, 2L, 2L)
id_col <- c(0L, 2L, 2L, 0L, 1L, 2L)
d <- seq(6)
d_mat <- reticulate::tuple(d, tuple(id_row, id_col, convert = F), covert = F)
m_csr <- scipy$sparse$csr_matrix(d_mat, shape = tuple(3L, 3L))
m_csc <- scipy$sparse$csc_matrix(d_mat, shape = tuple(3L, 3L))

r_from_csr <- py_to_r(m_csr)
## 3 x 3 sparse Matrix of class "dgRMatrix"
## [1,] 1 . 2
## [2,] . . 3
## [3,] 4 5 6

r_from_csc <- py_to_r(m_csc)
## 3 x 3 sparse Matrix of class "dgCMatrix"
## [1,] 1 . 2
## [2,] . . 3
## [3,] 4 5 6


# * load sample H3K27ac h5ad

# ** functions
#' @export
getAnnX <- function(ann, groupby, group) {
  allbarcodes <- py_to_r(ann$obs_names$values)
  allfeatures <- py_to_r(ann$var_names$values)
  groups <- py_to_r(ann$obs[[groupby]]$astype("str"))
  barcodes <- allbarcodes[groups == group]
  mat <- py_to_r(ann[barcodes, allfeatures]$X)
  rownames(mat) <- barcodes
  colnames(mat) <- allfeatures
  return(mat)
}

#' Generate Monocle3 object for cicero analysis
#'
#' In Monocle3, the matrix can be sparse or not.
#' But the row of the stored matrix will be genes or features
#' Also the matrix will be saved as dgCMatrix if sparse.
#' If dgRMatrix is used, will be transformed into dgCMatrix, which
#' needs to load Matrix package explicitly.
#'
#' Matrix will be binarzied if they have reads / counts or not,
#' and full zero cells and features will be filtered out.
c #'
#' FIXME:
#' > Warning messages:
#'  The `x` argument of `as_tibble.matrix()`
#'  must have unique column names if `.name_repair`
#'  is omitted as of tibble 2.0.0. Using compatibility `.name_repair`.
#'
#' @param mat matrix or sparse matrix with rownames (features) and
#'   colnames (cells)
#' @param cellinfo data.frame rownames is the same as the colnames of
#'   mat If NULL, will transform colnames of mat to data.frame.
#' @param peakinfo data.frame rownames is the same as the rownames of
#'   mat If NULL, will transform rownames of mat to data.frame.
#' @param isSparse bool for mat
#' @param isRowBarcode bool for mat
#' @param split characters for spliting peak string
#' @return monocle3::cell_data_set object
#' @export
getBinarizeMonocle3Obj <- function(mat,
                                   cellinfo = NULL,
                                   peakinfo = NULL,
                                   isSparse = T,
                                   isRowBarcode = T,
                                   split = ":|-") {
  bmat <- mat
  # binarization
  if (isSparse) {
    bmat@x[bmat@x > 0] <- 1
  } else {
    bmat[bmat > 0] <- 1
  }
  # features as rows
  if (isRowBarcode) {
    bmat <- t(bmat)
  }
  # remove full-zero rows and columns
  cols_zero <- colSums(bmat) == 0
  if (sum(cols_zero) > 0) {
    bmat <- bmat[, !cols_zero]
  }
  rows_zero <- rowSums(bmat) == 0
  if (sum(rows_zero) > 0) {
    bmat <- bmat[!rows_zero, ]
  }
  # set cellinfo and feainfo if they are NULL
  if (is.null(cellinfo)) {
    cellinfo <- data.frame(
      barcode = colnames(bmat),
      row.names = colnames(bmat)
    )
  }
  if (is.null(peakinfo)) {
    peakinfo <- strsplit(
      rownames(bmat),
      split = ":|-",
      perl = TRUE
    ) |>
      do.call(rbind, args = _) |>
      tibble::as_tibble() |>
      setNames(object = _, c("chr", "bp1", "bp2")) |>
      dplyr::mutate(
        bp1 = as.integer(bp1),
        bp2 = as.integer(bp2),
        site_peak = rownames(bmat)
      ) |>
      as.data.frame()
    rownames(peakinfo) <- peakinfo$site_peak
  }

  # convert dgR to dgC
  # need to load Matrix package
  if (class(bmat) == "dgRMatrix") {
    bmat <- as(bmat, "CsparseMatrix")
  }
  # to monocle3 CDS
  # repress warning below:
  # In monocle3::new_cell_data_set
  #    gene_metadata must contain a column verbatim named
  #    'gene_short_name' for certain functions.
  cds <- suppressWarnings(monocle3::new_cell_data_set(
    expression_data = bmat,
    cell_metadata = cellinfo,
    gene_metadata = peakinfo
  ))
  cds <- monocle3::detect_genes(cds)
  # Ideally we don't need this, but keep it to follow tutorial.
  cds <- cds[rowSums(exprs(cds)) != 0, ]
  return(cds)
}


# ** load whole matrix
ann <- ad$read_h5ad(
  file.path(projd, "12.DE", "out", "ann.H3K27ac.pmat.h5ad"),
  backed = "r"
)

# dgRMatrix
pmat <- py_to_r(ann$X)
# if no values, it does not conver to array
barcodes <- py_to_r(ann$obs_names$values)
features <- py_to_r(ann$var_names$values)
rownames(pmat) <- barcodes
colnames(pmat) <- features

# ** load partial matrix
# View of ann, and can be used to load the data fast.
groupby <- "annot.sc"
group <- "334 Microglia NN"

pmat_MGA <- getAnnX(ann, "annot.sc", "334 Microglia NN")
pmat_BAM <- getAnnX(ann, "annot.sc", "335 BAM NN")

cds_MGA <- getBinarizeMonocle3Obj(pmat_MGA)

# perform dimension reduction using monocle3
# Ideally we can use whatever method, but failed to use snapatac2 here.
# So keep the default one.
set.seed(2024)
cds_MGA <- monocle3::estimate_size_factors(cds_MGA) |>
  monocle3::preprocess_cds(
    cds = _,
    method = "LSI", num_dim = 50
  ) |>
  monocle3::reduce_dimension(
    cds = _, max_components = 2,
    reduction_method = "UMAP", preprocess_method = "LSI", cores = 1
  )
cds_cicero_MGA <- cicero::make_cicero_cds(
  cds = cds_MGA,
  reduced_coordinates = SingleCellExperiment::reducedDims(cds_MGA)$UMAP,
  k = 50, return_agg_info = T
)

# * run cicero by calling functions seperately
# load mm10 genome size
mouse_mm10_genome <- read.table(
  file = file.path(projd, "meta", "mm10.chrom.sizes.lite"),
  header = F, sep = "\t"
)

sample_genome <- subset(mouse_mm10_genome, V1 == "chr2")
sample_genome$V2[1] <- 10000000

# usually run with sample_num = 100
conns <- run_cicero(cds_cicero_MGA[[1]], sample_genome, sample_num = 10)
conns <- connes
# remove NA values
conns <- conns[!is.na(conns$coaccess), ]

# ** run BAM
cds_BAM <- getBinarizeMonocle3Obj(pmat_BAM) |>
  monocle3::estimate_size_factors(cds = _) |>
  monocle3::preprocess_cds(
    cds = _,
    method = "LSI", num_dim = 50
  ) |>
  monocle3::reduce_dimension(
    cds = _, max_components = 2,
    reduction_method = "UMAP", preprocess_method = "LSI", cores = 1
  )

cds_cicero_BAM <- cicero::make_cicero_cds(
  cds = cds_BAM,
  reduced_coordinates = SingleCellExperiment::reducedDims(cds_BAM)$UMAP,
  k = 50, return_agg_info = T
)

# After testing make_cicero_cds,
# let's run cicero for all the cells in one model.
# Under this way, we aim to find the enhancer-gene connections
# that show consistent regulations along all the cells / cell types.
# And this probably make correlation analysis
# using pseudobulk RNA-seq make more sense.
# But this would lose the
# cell-type specific cis-co-accessibility networks (CCANs)
# we can try the CCANs.

# * parallizaiton on cicero
# It seems that future will not copy all the data for independent sessions
# so this will save lots of memory
## for debug
generate_windows <- cicero:::generate_windows
get_genomic_range <- cicero:::get_genomic_range
calc_dist_matrix <- cicero:::calc_dist_matrix
find_distance_parameter <- cicero:::find_distance_parameter

cds <- cds_cicero_BAM[[1]]
window <- 5e+05
maxit <- 100
s <- 0.75
sample_num <- 10
distance_constraint <- 250000
distance_parameter_convergence <- 1e-22
max_elements <- 200
genomic_coords <- mouse_mm10_genome
genomic_coords <- genomic_coords[genomic_coords$V1 == "chr1", ]
genomic_coords$V2 <- 10000000
max_sample_windows <- 500

p_estimate_distance_parameter <- function(
    cds, window = 5e+05, maxit = 100, s = 0.75,
    sample_num = 100,
    distance_constraint = 250000,
    distance_parameter_convergence = 1e-22,
    max_elements = 200,
    genomic_coords = cicero::human.hg19.genome,
    max_sample_windows = 500) {
  grs <- generate_windows(window, genomic_coords)
  fData(cds)$chr <- gsub("chr", "", fData(cds)$chr)
  fData(cds)$bp1 <- as.numeric(as.character(fData(cds)$bp1))
  fData(cds)$bp2 <- as.numeric(as.character(fData(cds)$bp2))

  ## this part can be replaced by future apply
  distance_parameters <- future_lapply(seq_len(max_sample_windows), \(it) {
    win <- sample(seq_len(length(grs)), 1)
    GL <- "Error"
    win_range <- get_genomic_range(grs, cds, win)
    if (nrow(exprs(win_range)) <= 1) {
      return(NA)
    }
    if (nrow(exprs(win_range)) > max_elements) {
      return(NA)
    }
    dist_matrix <- calc_dist_matrix(win_range)
    distance_parameter <- find_distance_parameter(dist_matrix,
      win_range,
      maxit = maxit,
      null_rho = 0, s,
      distance_constraint = distance_constraint,
      distance_parameter_convergence = distance_parameter_convergence
    )
    if (!is(distance_parameter, "numeric")) {
      # which m
      return(NA)
    }
    return(distance_parameter)
  }, future.seed = T)

  ps <- unlist(distance_parameters)
  if (sum(is.na(ps)) == length(ps)) {
    stop("No distance_parameters calculated")
  }
  ps <- ps[!is.na(ps)]
  if (length(ps) < sample_num) {
    warning(paste0(
      "Could not calculate sample_num distance_parameters (",
      length(ps), " were calculated) - see ",
      "documentation details"
    ))
  }
  if (length(ps) > sample_num) {
    ps <- sample(ps, size = sample_num, replace = F)
  }
  return(ps)
}

## for debug
get_rho_mat <- cicero:::get_rho_mat
distance_parameter <- mean(ps)

p_generate_cicero_models <- function(
    cds, distance_parameter, s = 0.75, window = 5e+05,
    max_elements = 200,
    genomic_coords = cicero::human.hg19.genome) {
  grs <- generate_windows(window, genomic_coords)
  fData(cds)$chr <- gsub("chr", "", fData(cds)$chr)
  fData(cds)$bp1 <- as.numeric(as.character(fData(cds)$bp1))
  fData(cds)$bp2 <- as.numeric(as.character(fData(cds)$bp2))

  outlist <- future_lapply(seq_len(length(grs)),
    function(win) {
      GL <- "Error"
      win_range <- get_genomic_range(grs, cds, win)
      if (nrow(exprs(win_range)) <= 1) {
        return("Zero or one element in range")
      }
      if (nrow(exprs(win_range)) > max_elements) {
        return("Too many elements in range")
      }
      dist_matrix <- calc_dist_matrix(win_range)
      rho_mat <- get_rho_mat(
        dist_matrix, distance_parameter,
        s
      )
      vals <- as.matrix(exprs(win_range))
      cov_mat <- cov(t(vals))
      diag(cov_mat) <- diag(cov_mat) + 1e-04
      GL <- glasso::glasso(cov_mat, rho_mat)
      colnames(GL$w) <- row.names(GL$w) <- row.names(vals)
      colnames(GL$wi) <- row.names(GL$wi) <- row.names(vals)
      return(GL)
    },
    future.seed = T
  ) ## end of parallization

  names_df <- as.data.frame(grs)
  names(outlist) <- paste(names_df$seqnames, names_df$start,
    names_df$end,
    sep = "_"
  )
  outlist
}
