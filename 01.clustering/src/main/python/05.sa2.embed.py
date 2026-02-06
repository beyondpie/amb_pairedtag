import pyprojroot
import polars
import os
import sys
from typing import List, Optional
from logging import Logger

import numpy as np
import pandas as pd
from numba.core.errors import NumbaDeprecationWarning
from numba.core.errors import NumbaPendingDeprecationWarning
import warnings

warnings.simplefilter("ignore", category=NumbaDeprecationWarning)
warnings.simplefilter("ignore", category=NumbaPendingDeprecationWarning)

import snapatac2 as sa2
import scanpy as sc

# scanpy variable feature needs this
# import skmisc
# import importlib
# importlib.reload(sa2)


proj_dir = pyprojroot.here()
pydir = os.path.join(proj_dir, "package", "python")
sys.path.insert(0, pydir)
import myumap
from mylog import set_file_logger

# * configs
eps: float = 0.0001

logfnm: str = snakemake.log[0]
b2g_fnm: str = snakemake.input["b2g_fnm"]
allen_k8: str = snakemake.input["allen_k8"]
ann_h5ad: str = snakemake.input["ann_h5ad"]
out_ann: str = snakemake.output["ann"]
out_cellmeta: str = snakemake.output["cellmeta"]
npc: int = snakemake.params["npc"]
k_knn: int = snakemake.params["k_knn"]
knn_method: str = snakemake.params["knn_method"]
# use for umap
threads: int = snakemake.threads
cl: int = int(snakemake.wildcards["i"])
# clustering level
cll: str = snakemake.params["cll"]
use_k8: int = snakemake.params["use_k8"]
nfeature: int = snakemake.params["nfeature"]

# DEBUG
# work_dir = f"{proj_dir}/01.clustering"
# out_dir = f"{work_dir}/out/test"
# os.makedirs(out_dir, exist_ok=True)
# logfnm = f"{out_dir}/test.clustering.log"
# b2g_fnm = f"{work_dir}/src/test/resource/b2g_test.csv"
# ann_h5ad = os.path.join(work_dir, "src/test/resource",
#                         "scRNAseq_sa2_test.ann.h5ad")
# allen_k8 = f"{proj_dir}/meta/AIT21_k8_markers.txt"
# npc = 50
# k_knn = 50
# knn_method = "kdtree"
# prefix = f"RNA_k8_npc{npc}_k{k_knn}_L1"
# out_ann = os.path.join(out_dir, f"{prefix}_0.h5ad")
# out_cellmeta = os.path.join(out_dir, f"{prefix}_0.csv")
# cl: int = 1
# cll: str = "L2"
# use_k8: int = 0
# nfeature: int = 3000

if knn_method == "pynndescent":
    raise RuntimeError(f"SnapATAC2 with knn method {knn_method} has bug.")

# * set logger
logger: Logger = set_file_logger(logfnm, name="sa2_embed_with_scRNAseq")

logger.info(f"SnapATAC2 version: {sa2.__version__}.")
logger.info(f"Scanpy version: {sc.__version__}.")

# * load input
b2g: pd.DataFrame = pd.read_csv(b2g_fnm, sep=",", header=0)
with open(allen_k8, "r") as f:
    k8s = [g.strip() for g in f.readlines()]
ann: sa2.AnnData = sa2.read(filename=ann_h5ad, backed="r")

# * subset ann firstly
barcodes: List[str] = b2g[b2g[cll] == cl].barcode.tolist()
logger.info(f"Under clustering level {cll}, group {cl}: {len(barcodes)} cells.")

genes_ann: List[str] = ann.var_names

if use_k8 > 0:
    k8_ann: List[bool] = [True if g in k8s else False for g in genes_ann]
    logger.info(f"After intersect with k8_allen, {sum(k8_ann)} features left.")
    genes_ann = k8_ann
else:
    logger.info("Will select variable features instead of using k8.")

ann_var: sa2.AnnData = ann.subset(
    obs_indices=barcodes, var_indices=genes_ann, out=out_ann
)
ann.close()

# * select features with variations
ann_in_mem = ann_var.to_memory()
logger.info("Select features with variations.")
max_exp = ann_in_mem.X.max(axis=0).toarray()[0, :]
min_exp = ann_in_mem.X.min(axis=0).toarray()[0, :]
delta_exp = max_exp - min_exp
gene_with_var = delta_exp > eps

k8_ann: List[bool] = [
    True if g in k8s else False for g in ann_in_mem.var_names]
k8_selected = polars.Series(np.logical_and(gene_with_var, k8_ann))
logger.info(f"Default #k8 with variation is {k8_selected.sum()}.")

if use_k8 > 0:
    logger.info(f"Use k8 features as needed.")
    selected = k8_selected
else:
    logger.info(f"Select top {nfeature} variable features using Scanpy.")
    logger.info("Use Seurat_v3, which is default vst method in Seurat.")
    try:
        vg: Optional[pd.DataFrame] = sc.pp.highly_variable_genes(
            adata=ann_in_mem, n_top_genes=nfeature, flavor="seurat_v3",
            inplace=False)
        if vg is not None:
            selected = polars.Series(np.logical_and(gene_with_var, vg.highly_variable))
            logger.info(f"{sum(selected)} features are kept.")
        else:
            logger.error("Get None from scanpy highly_Variable_genes fn.")
            logger.info("Now use k8 genes as default feature for later analysis.")
            selected = k8_selected
    except ValueError as e:
        logger.error("Failed on sc.pp.highly_variable_genes. Now use k8 features as default.")
        logger.error(e)
        selected = k8_selected

ann_var.var["selected"] = selected
# * non-linear embed
logger.info("Run spectral embedding.")
sa2.tl.spectral(
    adata=ann_var, n_comps=npc, features="selected", weighted_by_sd=True,
    inplace=True
)
# * umap
logger.info("Run UMAP with given a and b.")
myumap.umap(
    adata=ann_var,
    n_comps=2,
    use_rep="X_spectral",
    key_added="umap",
    inplace=True,
    # use no seed for parallism
    # random_state=0,
    random_state=None,
    n_jobs=threads,
    a=1.8956,
    b=0.8006,
    init="spectral",
    metric="euclidean",
)
# * knn
logger.info(f"Run KNN with k={k_knn} and method={knn_method}.")
sa2.pp.knn(
    adata=ann_var,
    n_neighbors=k_knn,
    use_dims=None,
    use_rep="X_spectral",
    method=knn_method,
    inplace=True,
    random_state=0,
)

logger.info("Finish to save ann data with embedding, knn and umap.")
# * summarize results
barcodes = pd.DataFrame(
    data=ann_var.obs_names, dtype=str, index=ann_var.obs_names,
    columns=["barcode"]
)
umap = pd.DataFrame(
    data=ann_var.obsm["X_umap"],
    index=ann_var.obs_names,
    columns=["umap1", "umap2"],
    dtype=float,
)
cell_meta: pd.DataFrame = pd.concat([barcodes, umap], axis=1)
cell_meta.to_csv(out_cellmeta, sep=",", header=True, index=False)
logger.info("Finish to save cell_meta.")
logger.info("Close SnapATAC2 AnnData.")
ann_var.close()
