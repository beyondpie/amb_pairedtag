#' Calculate Overlap Score matrix for ATAC and RNA jointly embedding metaTable
#' @description
#' Based on the co-embeding clustering result, we summation the minmum scores of
#' the percentages of cells on either ATAC or RNA group in each cluster.
#' @param meta data.frame,
#' three columns as ident, atacCluster, rnaCluster defined as in the following
#' parameters
#' @param ident characters, name of column ident, i.e. the cluster Ids for
#' the co-embedded seurat, "coembed.idents" as default.
#' @param atacCol characters, name of column atac,
#' "MajorType"as default
#' @param rnaCol characters, name of column rna,
#' "ClusterName" as default
#' @return data.frame of numeric, atac by rna
#' @export
getOverlapMatrix <- function(meta,
                             ident = "coembed.idents",
                             atacCol = "MajorType",
                             rnaCol = "ClusterName") {
  ident2rna <- data.frame(idents = meta[[ident]], rna_label = meta[[rnaCol]])
  ident2rna <- ident2rna[stats::complete.cases(ident2rna), ]
  ident2atac <- data.frame(idents = meta[[ident]], atac_label = meta[[atacCol]])
  ident2atac <- ident2atac[stats::complete.cases(ident2atac), ]
  rnaTable <- table(ident2rna)
  atacTable <- table(ident2atac)
  rnaPct <- apply(rnaTable, 2, function(x) {
    x / sum(x)
  })
  atacPct <- apply(atacTable, 2, function(x) {
    x / sum(x)
  })
  rnaClusterName <- colnames(rnaPct)
  atacClusterName <- colnames(atacPct)
  calOvlpScoreElement <- function(t1l, t2l) {
    t1PctDF <- data.frame(rnaPct[, t1l])
    colnames(t1PctDF) <- "t1"
    t1PctDF$ident <- rownames(t1PctDF)
    t2PctDF <- data.frame(atacPct[, t2l])
    colnames(t2PctDF) <- "t2"
    t2PctDF$ident <- rownames(t2PctDF)
    comp <- plyr::join(t1PctDF, t2PctDF, by = "ident", type = "full")
    comp[is.na(comp)] <- 0
    comp$ident <- NULL
    comp <- t(comp)
    return(sum(apply(comp, 2, min)))
  }
  rna2atacType <- outer(rnaClusterName, atacClusterName,
    FUN = paste, sep = "|"
  )
  ovlpScore <- apply(rna2atacType, MARGIN = c(1, 2), FUN = function(i) {
    t <- unlist(strsplit(i, split = "|", fixed = TRUE))
    calOvlpScoreElement(t[1], t[2])
  })
  rownames(ovlpScore) <- rnaClusterName
  colnames(ovlpScore) <- atacClusterName
  return(t(ovlpScore))
}

#' mapGmatToGeneName
#' @return data.table, cell by gene
#' @export
mapGmatToGeneName <- function(gmat,
                              gencode,
                              cells,
                              fn = base::rowMeans) {
  geneIds <- colnames(gmat)
  geneNames <- gencode[geneIds, gene_name]
  dupIndex <- duplicated(geneNames)
  dupNames <- unique(geneNames[dupIndex])

  gdupl <- lapply(dupNames, function(g) {
    v <- gmat[, .(fn(do.call(cbind, .SD))),
      .SDcols = which(geneNames %in% g)
    ]
    colnames(v) <- g
    return(v)
  })
  gdup <- do.call(cbind, gdupl)

  gdup$cell <- cells
  data.table::setkey(gdup, cell)
  guniq <- gmat[, .SD, .SDcols = which(!(geneNames %in% dupNames))]
  colnames(guniq) <- geneNames[!(geneNames %in% dupNames)]
  guniq$cell <- cells
  data.table::setkey(guniq, cell)
  r <- guniq[gdup, on = .(cell = cell)]
  r[, cell := NULL]
  return(r)
}

#' @export
getAllenSubclassGeneMarkers <- function(allenMeta,
                                        all.genes,
                                        subclass = "LA-BLA-BMA-PA Glut") {
  index.allen <- allenMeta$upl2 %in% subclass
  # gene from name
  name.pool <- as.vector(
    stringr::str_split(subclass, pattern = " ", simplify = TRUE)
  )
  name.index <- name.pool %in% all.genes
  name.markers <- if (sum(name.index) > 0) {
    message(subclass, " has match genes in its name.")
    name.pool[name.index]
  } else {
    message(subclass, " has no matched gene in its name.")
    NULL
  }

  # top.markers
  top.markers.pool <- as.vector(stringr::str_split(
    allenMeta$top.markers[index.allen],
    pattern = ",",
    simplify = TRUE
  ))
  top.markers.pool <- top.markers.pool[nchar(top.markers.pool) > 0]
  top.markers.pool <- top.markers.pool[top.markers.pool %in% all.genes]
  top.markers <- names(sort(table(top.markers.pool), decreasing = TRUE))
  message(subclass, " has ", length(top.markers), " top markers.")

  # combo markers
  combo.markers.pool <- as.vector(stringr::str_split(
    allenMeta$combo.markers[index.allen],
    pattern = ",",
    simplify = TRUE
  ))
  combo.markers.pool <- combo.markers.pool[
    nchar(combo.markers.pool) > 0
  ]
  combo.markers.pool <- combo.markers.pool[
    combo.markers.pool %in% all.genes
  ]
  combo.markers <- names(sort(table(combo.markers.pool),
    decreasing = TRUE
  ))
  message(subclass, " has ", length(combo.markers), " combo markers.")
  # combine all the signals
  signal.nms <- unique(c(name.markers, top.markers, combo.markers))
  return(signal.nms)
}

#' @export
getIntOvlpMat <- function(class = "nn",
                          k.anchor = 5,
                          k.feature = 200,
                          intEmbed = "cca",
                          intRef = "notatac",
                          featuresFrom = "notatac",
                          atacEmbed = "pca",
                          annot1 = "L3",
                          annot2 = "l2",
                          r = 0.8) {
  rdir <- file.path(
    here::here(), "03.integrate",
    paste("int", class, "L3", sep = "."), "out"
  )
  fnm <- stringr::str_glue(
    "kAnchor{k.anchor}-kFeature{k.feature}",
    "{intEmbed}_Int-ref_{intRef}",
    "feat_{featuresFrom}",
    "atacEmbed_{atacEmbed}",
    "annot-{annot1}_{annot2}.rds",
    .sep = "-"
  )
  # load data
  result <- readRDS(file.path(rdir, fnm))
  mat <- result$intgnClusterList[[paste0("r", r)]][["ovlpMat"]]
  return(mat)
}

# ===========================================
# Add Seurat v5 integration related functions
# ===========================================
#' @export
convertAnn2Seurat5 <- function(annfnm,
                               modality,
                               group = "X",
                               outdir,
                               overwrite = TRUE,
                               assay = "RNA",
                               isLogNorm = TRUE,
                               removeCounts = FALSE,
                               addObsMeta = FALSE,
                               usePython = "") {
  if (dir.exists(outdir) && (!overwrite)) {
    message(outdir, " exist and no overwrite.")
  } else {
    message("read anndata from hdf5: ", annfnm)
    xlognorm <- BPCells::open_matrix_anndata_hdf5(
      path = annfnm, group = group
    )
    message("write matrix to dir: ", outdir)
    BPCells::write_matrix_dir(mat = xlognorm, dir = outdir, overwrite = overwrite)
  }
  d <- BPCells::open_matrix_dir(outdir)
  s5 <- Seurat::CreateSeuratObject(counts = d, assay = assay)
  if (isLogNorm) {
    s5 <- Seurat::SetAssayData(
      object = s5, slot = "data", new.data = d
    )
    if (removeCounts) {
      message("remove counts slot in assay: ", assay)
      message("This may affect FindVariableFeatures function.")
      s5[[assay]]$counts <- NULL
    }
  }
  s5$modality <- modality
  if (addObsMeta) {
    message("Add obs_meta to Seurat.")
    if (!grepl("python", usePython)) {
      warning(usePython, ": no python detected.")
    } else {
      reticulate::use_python(usePython)
      ad <- reticulate::import("anndata", convert = FALSE)
      allen_ann <- ad$read_h5ad(filename = annfnm, backed = "r")
      obs_meta <- allen_ann$obs
      obsMeta <- reticulate::py_to_r(obs_meta)
      attr(obsMeta, "pandas.index") <- NULL
      if ("cl" %in% colnames(obsMeta)) {
        cls <- obsMeta$cl
        if (is.factor(cls)) {
          message("Hard code 'cl' column from factor to int.")
          clInt <- as.integer(levels(cls)[cls])
          obsMeta$cl <- clInt
        }
      } # end of transform cl to Int
      obsMeta <- obsMeta[colnames(s5), ]
      s5 <- Seurat::AddMetaData(s5, metadata = obsMeta)
    } # end of handling of adding obs meta in details
  } # end of adding obs meta
  return(s5)
}

#' depend on get.downsample.fun in the same file
#' @export
downsampleSeurat <- function(seu,
                             groupBy = "cluster_id",
                             minNum = 1000,
                             maxNum = 2000) {
  fn.dp <- get.downsample.fun(
    minNum = minNum, maxNum = maxNum
  )
  allcells <- colnames(seu)
  ## NOTE: seu[[groupBy]] will return data.frame with one column
  labels <- seu@meta.data[[groupBy]]
  if (is.null(labels)) {
    warning(groupBy, " is not int meta.data.")
    return(seu)
  }
  dp.cells <- fn.dp(index = allcells, labels = labels)
  subset(seu, cells = dp.cells)
}

#' @export
toSeuratInMemory <- function(seu, slot = "data", assay = "RNA",
                             removeCounts = FALSE) {
  mat <- as(seu[[assay]][[slot]], Class = "dgCMatrix")
  meta <- seu@meta.data
  message("transform to in-memory seurat object with slot: ", slot)
  r <- Seurat::CreateSeuratObject(
    counts = mat, assay = assay,
    meta.data = meta
  )
  if (slot == "data") {
    r <- Seurat::SetAssayData(
      object = r,
      slot = "data", new.data = mat
    )
    if (removeCounts) {
      message("remove counts slot for in-memory seurat.")
      r[[assay]]$counts <- NULL
    }
  }
  return(r)
}

#' @export
isOnDiskMat.Seurat <- function(seu, onlayer = "counts") {
  t <- SeuratObject::LayerData(seu, layer = onlayer)
  if (attr(class(t), "package") == "BPCells") {
    TRUE
  } else {
    FALSE
  }
}

#' depends on isOnDiskMat.Seurat
#' @export
calVarOfFea.Seurat <- function(seu, onlayer = "counts") {
  mat <- SeuratObject::LayerData(seu, layer = onlayer)
  if (isOnDiskMat.Seurat(seu, onlayer)) {
    tmp <- BPCells::matrix_stats(
      matrix = mat, row_stats = "variance", col_stats = "none"
    )
    tmp$row_stats["variance", ]
  } else {
    sparseMatrixStats::rowVars(mat) |>
      x => stats::setNames(x, rownames(seu))
  }
}

#' @export
getvf.Seurat <- function(seu,
                         onlayer = "count",
                         eps = 0.0001) {
  nall <- nrow(seu)
  x <- calVarOfFea.Seurat(seu, onlayer)
  vf <- names(x)[x > eps]
  nf <- length(vf)
  message(
    "From ", nall,
    " features, detect ", nf, " variational features."
  )
  return(vf)
}

#' depends on isOnDiskMat.Seurat, calVarOfFea.Seuat
#' @export
setVariableFeatures <- function(ref,
                                query,
                                features = NULL,
                                onlayer = "counts",
                                nfeatures = 2000,
                                on = c("ref", "all"),
                                eps = 0.0001) {
  message("Get variational features of ref.")
  vfref <- getvf.Seurat(ref, onlayer, eps)
  message("Get variational features of query.")
  vfquery <- getvf.Seurat(query, onlayer, eps)


  vf <- if (!is.null(features)) {
    nf <- length(features)
    message("Detect ", nf, " features.")
    message("Ignore variablefeatures then.")
    intersect(features, intersect(vfref, vfquery))
  } else {
    useOn <- match.arg(on)
    message("Set variable features based on: ", useOn)

    if (useOn == "ref") {
      refVF <- if (length(x = Seurat::VariableFeatures(ref)) == 0) {
        message("No variable features in ref. Calculating...")
        Seurat::FindVariableFeatures(ref, nfeatures = nfeatures) |>
          x => Seurat::VariableFeatures(x)
      } else {
        Seurat::VariableFeatures(ref)
      }
      intersect(refVF, intersect(vfref, vfquery))
    } else {
      Seurat::SelectIntegrationFeatures(
        object.list = list(ref, query),
        nfeatures = nfeatures,
        fvf.nfeatures = nfeatures,
        verbose = TRUE
      )
    }
  }
  message("Get ", length(vf), " features finally.")
  return(vf)
}

#' @export
getPredictLabel <- function(tfquery) {
  labels <- tfquery$predicted.id
  scores <- tfquery$predicted.id.score
  cells <- colnames(tfquery)
  r <- data.frame(
    barcode = cells,
    predict = labels,
    score = scores
  )
  rownames(r) <- cells
  return(r)
}

#' @export
getPredictScoreMat <- function(tfquery) {
  predictAssay <- tfquery@assays$prediction.score.id
  score <- predictAssay@data
  return(score)
}

#' getMetaCol.Seurat
#' 
#' use repeat colnm to get array
#' instead of data.frame with one column
#' @return array with name
#' @export
getMetaCol.Seurat <- function(seu, colnm) {
  g <- seu[[colnm]][[colnm]]
  names(g) <- colnames(seu)
  return(g)
}


#' tfmat: transfer label score: ref by query
#' @export
prepareDotPlot4TransferLabel <- function(tfmat,
                                         refOrder = NULL,
                                         names = c("row", "column", "score"),
                                         ignoreEmptyRef = TRUE) {
  maxScore <- getMaxColScoreWithName(mat = tfmat)
  query2ref <- data.frame(
    query = colnames(tfmat),
    ref = names(maxScore),
    row.names = colnames(tfmat)
  )
  if (ignoreEmptyRef) {
    message("remove refs not having query mapped to.")
    tfmat <- tfmat[rownames(tfmat) %in% query2ref$ref, ]
  }
  if (is.null(refOrder)) {
    message("refOrder is null, will use default numeric order for it.")
    refOrder <- rownames(tfmat)
    refOrder <- refOrder[order(as.integer(refOrder))]
  } else {
    refOrder <- refOrder[refOrder %in% rownames(tfmat)]
  }
  queryOrder <- query2ref$query[
    order(factor(query2ref$ref, levels = refOrder))
  ]

  meltmat <- to3.matrix(tfmat, names)
  meltmat[, 1] <- factor(meltmat[, 1], levels = refOrder)
  meltmat[, 2] <- factor(meltmat[, 2], levels = queryOrder)
  # reduce size of meltmat
  meltmat <- meltmat[meltmat[, 3] > 0, ]
  return(meltmat)
}

#' @export
loadAllenAnnot <- function(fnm) {
  r <- data.table::fread(
    file = fnm, sep = "\t",
    header = TRUE, data.table = FALSE
  )
  rownames(r) <- r$cl
  return(r)
}
