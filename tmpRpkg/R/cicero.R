#' Get the full sparse matrix from AnnData
#'
#' @return sparse matrix (dgR / dgC depends on AnnData)
#' with barcodes as rownames and features as colnames
#' @export
getFullAnnX <- function(ann) {
  allbarcodes <- reticulate::py_to_r(ann$obs_names$values)
  allfeatures <- reticulate::py_to_r(ann$var_names$values)
  # ann$X will not be transformed into R object
  mat <- reticulate::py_to_r(ann[allbarcodes, allfeatures]$X)
  rownames(mat) <- allbarcodes
  colnames(mat) <- allfeatures
  return(mat)
}

#' Get part of sparse matrix from AnnData
#' @return sparse matrix (dgR / dgC depends on AnnData)
#' with barcodes as rownames and features as colnames
#' @export
getPartAnnX <- function(ann, groupby, group) {
  allbarcodes <- reticulate::py_to_r(ann$obs_names$values)
  allfeatures <- reticulate::py_to_r(ann$var_names$values)
  groups <- reticulate::py_to_r(ann$obs[[groupby]]$astype("str"))
  barcodes <- allbarcodes[groups == group]
  mat <- reticulate::py_to_r(ann[barcodes, allfeatures]$X)
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
#' Matrix will be binarzied if they have reads orcounts or not,
#' and full zero cells and features will be filtered out.
#'
#' FIXME:
#' Warning messages:
#'  The x argument of as_tibble.matrix()
#'  must have unique column names if .name_repair
#'  is omitted as of tibble 2.0.0. Using compatibility .name_repair.
#'
#' @param mat matrix or sparse matrix with rownames (features) and
#'   colnames (cells)
#' @param cellinfo data.frame rownames is the same as the colnames of
#'   mat If NULL, will transform colnames of mat to data.frame.
#' @param peakinfo data.frame rownames is the same as the rownames of
#'   mat If NULL, will transform rownames of mat to data.frame.
#' @param isSparse bool for mat
#' @param isRowBarcode bool for mat
#' @param splitby characters for spliting peak string
#' @return monocle3::cell_data_set object
#' @export
getBinarizeMonocle3Obj <- function(mat,
                                   cellinfo = NULL,
                                   peakinfo = NULL,
                                   isSparse = T,
                                   isRowBarcode = T,
                                   splitby = ":|-") {
  bmat <- mat
  # binarization
  if (isSparse) {
    bmat@x[bmat@x > 0] <- 1
  } else {
    bmat[bmat > 0] <- 1
  }
  # features as rows
  if (isRowBarcode) {
    bmat <- Matrix::t(bmat)
  }
  # remove full-zero rows and columns
  cols_zero <- Matrix::colSums(bmat) == 0
  if (sum(cols_zero) > 0) {
    bmat <- bmat[, !cols_zero]
  }
  rows_zero <- Matrix::rowSums(bmat) == 0
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
      split = splitby,
      perl = TRUE
    ) |>
      do.call(rbind, args = _) |>
      tibble::as_tibble(x = _) |>
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
  cds <- cds[Matrix::rowSums(exprs(cds)) != 0, ]
  return(cds)
}

#' @export
generate_windows <- cicero:::generate_windows

#' @export
get_genomic_range <- cicero:::get_genomic_range

#' @export
calc_dist_matrix <- cicero:::calc_dist_matrix

#' @export
find_distance_parameter <- cicero:::find_distance_parameter

#' @export
get_rho_mat <- cicero:::get_rho_mat

#' @export
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

#' @export
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

#' Read positive proximal-distal connections of bedpe format.
#' @param fnm \link{ciceroPosProximalDistalConnfnm}
#' @return data.frame with columns defined
#' @export
loadppdc <- function(fnm) {
  data.table::fread(
    file = fnm, sep = "\t", header = F,
    data.table = F
  ) |>
    setNames(object = _, nm = c(
      "pchr", "pstart", "pend", "dchr", "dstart", "dend",
      "gene|enhancer", "coacc", "n1", "n2"
    ))
}

#' Read all the subclasses' positive proximal-enhancer connections
#'
#' @param ppecd \link{snATACppecd}
#' @return list of data.frame
#' @export
loadAllppec <- function(ppecd) {
  scs <- list.files(ppecd, full.names = F, no.. = T, include.dirs = F) |>
    x => basename(x) |>
    x => strsplit(x, split = ".", fixed = T) |>
    x => vapply(x, \(i) i[1], "a") |>
    unique()
  lapply(scs, \(sc) {
    loadppdc(fnm = str_glue("{ppecd}/{sc}.ppec.bedpe"))
  }) |> setNames(object = _, nm = scs)
}

#' Get distal2gene (by most correlations) from correlation of
#' pseudobulk-level RNA-seq and ATAC-seq.
#'
#' @param fnm \link{ppdcWithCorfnm}
#' @param genes vector of genes to keep, can be NULL
#' @param distals vector of CREs to keep, can be NULL
#' @return data.frame
#' two columns: distal and gene
#' @export
loadCRE2Gene <- function(fnm, genes = NULL, distals = NULL) {
  ppdcWithCor <- data.table::fread(
    file = fnm,
    header = T, sep = "\t", data.table = F
  )
  genes <- vapply(ppdcWithCor$conns, function(i) {
    unlist(strsplit(i, split = "|", fixed = TRUE))[1]
  }, FUN.VALUE = "a")
  distals <- vapply(ppdcWithCor$conns, function(i) {
    unlist(strsplit(i, split = "|", fixed = TRUE))[2]
  }, FUN.VALUE = "peak1")
  ppdcWithCor$gene <- genes
  ppdcWithCor$distal <- distals
  if (!is.null(genes)) {
    ppdcWithCor <- subset(ppdcWithCor, gene %in% genes)
  }
  if (!is.null(distals)) {
    ppdcWithCor <- subset(ppdcWithCor, distal %in% distals)
  }
  ppdcWithCor |>
    group_by(distal) |>
    summarise(geneByPeak = gene[which.max(pcc)]) |>
    as.data.frame() |>
    x => `rownames<-`(x, x$distal)
}
