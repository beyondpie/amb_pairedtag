library(BPCells)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(stringr)
library(purrr)
library(reticulate)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
reticulate::use_condaenv(
  condaenv = "sa2",
  conda = "/home/szu/miniforge3/bin/mamba",
  required = TRUE
)
ad <- reticulate::import(module = "anndata")

projdir <- here::here()
workdir <- file.path(projdir, "03.integration")
outdir <- file.path(projdir, "data/allen_seurat5")

transformAnnData2Seurat5 <- function(annfnm, outs5fnm,
                                     group = "X") {
  if (!file.exists(annfnm)) {
    stop(paste(annfnm, "does not exist."))
  }
  outs5dir <- dirname(outs5fnm)
  if (!dir.exists(outs5dir)) {
    message("Create Seurat5 dir: ", outs5dir)
    dir.create(path = outs5dir)
  }
  if (file.exists(outs5fnm)) {
    message(outs5fnm, " exist and will delete it first.")
    file.remove(outs5fnm)
  }
  message("Load matrix from the group: ", group)
  mat <- BPCells::open_matrix_anndata_hdf5(
    path = annfnm,
    group = group
  ) |> BPCells::convert_matrix_type(
    matrix = _, type = "uint32_t"
  )
  message(mat@dim[1], " genes as rows.")
  message(mat@dim[2], " cells as columns.")
  outmatdir <- file.path(outs5dir, "_mat")
  message("write matrix to ", outmatdir)
  if (dir.exists(outmatdir)) {
    message(outmatdir, " exists and remove it.")
    unlink(outmatdir)
  }
  BPCells::write_matrix_dir(
    mat = mat, dir = outmatdir,
    overwrite = TRUE
  )
  message("Load cell meta data.")
  ann <- ad$read_h5ad(filename = annfnm, backed = "r")
  cellmeta <- ann$obs
  d <- BPCells::open_matrix_dir(outmatdir)
  s5 <- Seurat::CreateSeuratObject(
    counts = d, assay = "RNA",
    meta.data = cellmeta
  )
  message("save seurat to rds file.")
  saveRDS(s5, file = outs5fnm)
  return(s5)
}

# * read Paired-Tag to Allen Region
p2r <- data.table::fread(
  file = file.path(workdir, "src/main/resource",
    "Allen_reference_data_locator_ZW_20240131.csv"),
  header = TRUE,
  sep = ",",
  data.table = FALSE
)

pairedtagRegions <- p2r[,1]
allen_10xv2_dir <- gsub("Allen 10x v2 \\(", "", colnames(p2r)[4]) |>
  gsub("\\)", "", x = _)
allen_10xv3_dir <- gsub("Allen 10x v3 \\(", "", colnames(p2r)[5]) |>
  gsub("\\)", "", x = _)

getAllenfnmOfpt <- function(rowid) {
  r <- p2r[rowid, 1]
  a1 <- if(is.na(p2r[rowid, 4])) {
    NA
  } else {
    file.path(allen_10xv2_dir, gsub(" ", "", str_split_1(p2r[rowid, 4], ";")))
  }
  a2 <- if(is.na(p2r[rowid, 5])) {
    NA
  } else {
    file.path(allen_10xv3_dir, gsub(" ", "", str_split_1(p2r[rowid, 5], ";")))
  }
  a <- c(a1, a2) |> x => x[!is.na(x)]
  data.frame(ptr = r, path = a)
}

ptr2path <- do.call(rbind, lapply(seq_len(nrow(p2r)), getAllenfnmOfpt))
all(vapply(ptr2path[,2], file.exists, TRUE))

outs5fnms <- vapply(seq_len(nrow(ptr2path)), function(i) {
  annpath <- ptr2path[i, 2]
  prefix <- gsub(".h5ad", "", basename(annpath))
  file.path(outdir, prefix, paste0(prefix, ".rds"))
}, "apath")

ptr2path["seurats5"] <- outs5fnms
write.table(ptr2path,
  file = file.path(workdir,
    "src/main/resource", "pairedtagRegion2AllenFile.csv"),
  sep = ",",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE)

# transform each file in ptr2path to Seurat object.
r <- lapply(seq_len(nrow(ptr2path)), function(i) {
  annpath <- ptr2path[i, 2]
  prefix <- gsub(".h5ad", "", basename(annpath))
  outs5fnm <- file.path(outdir, prefix, paste0(prefix, ".rds"))
  message("from ", annpath)
  message("to ", outs5fnm)
  transformAnnData2Seurat5(annfnm = annpath,
    bouts5fnm = outs5fnm, group = "layers/rawcount")
})

# * get Seurat by merging all the regions
pt2allenr <- data.table::fread(file = pt2allenRegionfnm,
  header = TRUE, sep = ",", data.table = FALSE)
all10xv3fnms <- with(pt2allenr, seurats5[grepl("10xv3", seurats5)])
names(all10xv3fnms) <- gsub(".rds", "", basename(all10xv3fnms))
all10xv3S5s <- lapply(all10xv3fnms, readRDS)
names(all10xv3S5s) <- names(all10xv3fnms)
ptregions <- with(pt2allenr, ptr[grepl("10xv3", seurats5)])
ncells <- vapply(all10xv3S5s, ncol, 1)

s510xv3 <- scCustomize::Merge_Seurat_List(
  list_seurat = all10xv3S5s)
s510xv3$barcode <- gsub("^_", "", colnames(s510xv3))

# downsample by region
ndp <- 250000
nr <- length(unique(ptregions))
ndpavg <- ceiling(ndp / nr)
barcodes <- s510xv3$barcode
regions <- lapply(seq_along(all10xv3S5s), \(i) {
  rep(ptregions[i], ncells[i]) }) |>
  do.call(what  = "c", args = _)
allenb2r <- data.frame(
  barcode = barcodes,
  region = regions
)

sub_allenb2r <- allenb2r |>
  dplyr::group_by(region) |>
  dplyr::slice_sample(n = ndpavg)

sub_s510xv3 <- subset(s510xv3,
  subset = barcode %in% sub_allenb2r$barcode)

sub_s510xv3 <- JoinLayers(sub_s510xv3)
# about 5 minutes, 30G
sub_s510xv3[["RNA"]]$counts <- as(
  object = sub_s510xv3[["RNA"]]$counts, Class = "dgCMatrix")
sub_s510xv3 <- SeuratObject::RenameCells(
  sub_s510xv3, new.names = sub_s510xv3$barcode)
# takes about 10 minutes with 3G on disk.
# maybe we can save the JoinLayers
saveRDS(sub_s510xv3, suballens5fnm)

# downsample on cl
ndp <- 100000

