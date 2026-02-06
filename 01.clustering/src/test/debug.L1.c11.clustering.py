
import sys; sys.path.extend(['/home/szu/miniforge3/envs/sa2/lib/python3.11/site-packages', '/home/szu/.cache/snakemake/snakemake/source-cache/runtime-cache/tmpxmh3cjlc/file/mnt/tscc2/szu/projects/amb_pairedtag/01.clustering/src/main/python', '/mnt/tscc2/szu/projects/amb_pairedtag/01.clustering/src/main/python']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95W\t\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94(\x8cO/projects/ps-renlab2/szu/projects/amb_pairedtag/meta/pt.barcode.meta.withL1.csv\x94\x8c\\/projects/ps-renlab2/szu/projects/amb_pairedtag/01.clustering/out/scRNAseq_sa2_afqc.ann.h5ad\x94\x8cI/projects/ps-renlab2/szu/projects/amb_pairedtag/meta/AIT21_k8_markers.txt\x94e}\x94(\x8c\x06_names\x94}\x94(\x8c\x07b2g_fnm\x94K\x00N\x86\x94\x8c\x08ann_h5ad\x94K\x01N\x86\x94\x8c\x08allen_k8\x94K\x02N\x86\x94u\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x18\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x1e)}\x94\x8c\x05_name\x94h\x18sNt\x94bh\x19h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x19sNt\x94bh\x10h\nh\x12h\x0bh\x14h\x0cub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94(\x8ce/mnt/tscc2/szu/projects/amb_pairedtag/01.clustering/out/afqc/L2/anns/RNA_vg3000_npc30_k50_L1_c11.h5ad\x94\x8ch/mnt/tscc2/szu/projects/amb_pairedtag/01.clustering/out/afqc/L2/RNA_vg3000_npc30_k50_L1_c11.cellmeta.csv\x94\x8ce/mnt/tscc2/szu/projects/amb_pairedtag/01.clustering/out/afqc/L2/flag/RNA_vg3000_npc30_k50_L1_c11.done\x94e}\x94(h\x0e}\x94(\x8c\x03ann\x94K\x00N\x86\x94\x8c\x08cellmeta\x94K\x01N\x86\x94\x8c\x04flag\x94K\x02N\x86\x94uh\x16]\x94(h\x18h\x19eh\x18h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x18sNt\x94bh\x19h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x19sNt\x94bh1h,h3h-h5h.ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94(K\x1eK2\x8c\x06kdtree\x94\x8c\x02L1\x94K\x00M\xb8\x0be}\x94(h\x0e}\x94(\x8c\x03npc\x94K\x00N\x86\x94\x8c\x05k_knn\x94K\x01N\x86\x94\x8c\nknn_method\x94K\x02N\x86\x94\x8c\x03cll\x94K\x03N\x86\x94\x8c\x06use_k8\x94K\x04N\x86\x94\x8c\x08nfeature\x94K\x05N\x86\x94uh\x16]\x94(h\x18h\x19eh\x18h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x18sNt\x94bh\x19h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x19sNt\x94bhHK\x1ehJK2hLhDhNhEhPK\x00hRM\xb8\x0bub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94\x8c\x0211\x94a}\x94(h\x0e}\x94\x8c\x01i\x94K\x00N\x86\x94sh\x16]\x94(h\x18h\x19eh\x18h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x18sNt\x94bh\x19h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x19sNt\x94b\x8c\x01i\x94haub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c\x04/tmp\x94e}\x94(h\x0e}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94uh\x16]\x94(h\x18h\x19eh\x18h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x18sNt\x94bh\x19h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x19sNt\x94bhxK\x01hzK\x01h|huub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94\x8cc/mnt/tscc2/szu/projects/amb_pairedtag/01.clustering/out/afqc/L2/log/RNA_vg3000_npc30_k50_L1_c11.log\x94a}\x94(h\x0e}\x94h\x16]\x94(h\x18h\x19eh\x18h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x18sNt\x94bh\x19h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x19sNt\x94bub\x8c\x06config\x94}\x94(\x8c\x04b2gf\x94\x8cO/projects/ps-renlab2/szu/projects/amb_pairedtag/meta/pt.barcode.meta.withL1.csv\x94\x8c\x07annh5ad\x94\x8c\\/projects/ps-renlab2/szu/projects/amb_pairedtag/01.clustering/out/scRNAseq_sa2_afqc.ann.h5ad\x94\x8c\x07allenk8\x94\x8cI/projects/ps-renlab2/szu/projects/amb_pairedtag/meta/AIT21_k8_markers.txt\x94\x8c\x06outdir\x94\x8c\x0bout/afqc/L2\x94\x8c\x03npc\x94K\x1e\x8c\x05k_knn\x94K2\x8c\x04ncpu\x94K\x01\x8c\x10clustering_level\x94hE\x8c\x06use_k8\x94K\x00\x8c\x08nfeature\x94M\xb8\x0b\x8c\x07nleiden\x94K\n\x8c\x06nsilht\x94K\x06\x8c\x07r_start\x94G?\xb9\x99\x99\x99\x99\x99\x9a\x8c\x05r_end\x94K\x02\x8c\rnsample_silht\x94MP\xc3\x8c\x0cnsample_umap\x94MP\xc3\x8c\x0cconda_prefix\x94\x8c\x14/home/szu/miniforge3\x94u\x8c\x04rule\x94\x8c\tsa2_embed\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8cC/mnt/tscc2/szu/projects/amb_pairedtag/01.clustering/src/main/python\x94ub.');

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
k8_ann: List[str] = [True if g in k8s else False for g in genes_ann]

if use_k8 > 0:
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
