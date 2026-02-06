#' @export
setHeatmap <- function(mat,
                       ha_row = NULL,
                       ha_col = NULL,
                       low.val.col = quantile(mat, 0.05, na.rm = TRUE),
                       high.val.col = quantile(mat, 0.95, na.rm = TRUE),
                       col_fun = circlize::colorRamp2(
                         seq(low.val.col, high.val.col, length = 60),
                         viridis::viridis(60)
                       ),
                       showRowNames = TRUE,
                       showColNames = FALSE,
                       fontsize = 6,
                       highqtl = 0.99,
                       cluster_columns = FALSE,
                       cluster_rows = FALSE,
                       use_raster = TRUE,
                       legend_show = TRUE,
                       legend_title = latex2exp::TeX(r"($\log$(CPM+1))"),
                       legend_direct = "horizontal") {
  legend_labels = c(
    round(low.val.col, 1),
    round(high.val.col, 1))
  p <- ComplexHeatmap::Heatmap(
    matrix = mat,
    col = col_fun,
    cluster_columns = cluster_columns,
    cluster_rows = cluster_rows,
    show_row_names = showRowNames,
    row_names_gp = grid::gpar(fontsize = fontsize),
    show_column_names = showColNames,
    top_annotation = ha_col,
    left_annotation = ha_row,
    use_raster = use_raster,
    show_heatmap_legend = legend_show,
    heatmap_legend_param = list(
      title = legend_title,
      at = c(low.val.col, high.val.col),
      labels = legend_labels,
      direction = legend_direct
    )
  )
  return(p)
}

#' Combine two ComplexHeatmaps.
#' @param chm1 ComplextHeatmap 1
#' @param chm2 ComplextHeatmap 2
#' @param legend_direct string, "bottom" by default
#' @param legend_merge boolean, TRUE by default
#' @return ComplextHeatmap
#' @export
combineHeatmap <- function(chm1, chm2,
                           legend_direct = "bottom",
                           legend_merge = TRUE) {
  ComplexHeatmap::draw(
    chm1 + chm2,
    heatmap_legend_side = legend_direct,
    merge_legend = legend_merge)
}

#' Group columns by hierarchical clustering
#' 
#' @param mat A matrix to be clustered
#' matrix should be in the form of samples x features
#' @param dmethod A method to calculate distance
#' @param hmethod A method to calculate hierarchical clustering
#' @param k A number of clusters
#' @param sumfn A function to summarize the grouped columns
#' in row-wise (will transform the mat during innter calculation). 
#' Default is rowMeans
#' @return list of enriched matrix, group (vector), and hclust object
#' @export
groupColsByHierarchical <- \(mat, dmethod = "pearson", hmethod = "ward.D2", k = 50, sumfn = mean) {
  res.cor <- stats::cor(mat, method = dmethod)
  res.dist <- stats::as.dist(1 - res.cor)
  hclustmat <- stats::hclust(res.dist, method = hmethod)
  g <- stats::cutree(hclustmat, k = k) |>
    x => paste0("g", x)
  enrich <- t(mat) |> as.data.frame()  |>
    mutate(group = g) |> 
    group_by(group) |>
    summarize(across(where(is.numeric), sumfn), 
              .groups = "drop") |>
              as.data.frame()
  rownames(enrich) <- enrich$group
  enrich <- enrich |> select(-group)
  return(list(enrich = enrich, group = g, hclust = hclustmat))
}

#' Get order of rows by orders of columns
#' @param enrich data.frame
#' @param topk topk to sort the rows
#' @return vector order of rows
#' @export
getOrderOfRowsByOrdersOfCols <- \(enrich, topk = 1) {
  topCols <- vapply(seq_len(nrow(enrich)), \(i) {
    order(unlist(enrich[i, ]), decreasing = T)[1:topk]
  }, FUN.VALUE = rep(1, topk))
  if (topk == 1) {
    topCols <- matrix(topCols, nrow = 1)
  }
  colnames(topCols) <- rownames(enrich)
  colMins <- colFn(min)
  avgRank <- colMins(topCols)
  rownames(enrich)[order(avgRank, decreasing = F)]
}

#' Get order of elements by group order
#' 
#' @param g vector of group info
#' @param group ordered unique names of g
#' @return vector order of g
#' @export
getOrderOfElementsByGroupOrder <- \(g, group) {
  r <- lapply(group, \(i) {
    seq_along(g)[which(g == i)]
  }) |> unlist()
  return(r)
}
