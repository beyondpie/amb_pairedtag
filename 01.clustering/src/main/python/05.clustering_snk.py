import pyprojroot
from sklearn.metrics import silhouette_score
import polars
import os
import sys
from typing import List, Tuple
from logging import Logger

import numpy as np
import pandas as pd
from numba.core.errors import NumbaDeprecationWarning
from numba.core.errors import NumbaPendingDeprecationWarning
import warnings
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)


import snapatac2 as sa2
# import importlib
# importlib.reload(sa2)


proj_dir = pyprojroot.here()
pydir = os.path.join(proj_dir, "package", "python")
sys.path.insert(0, pydir)
import myumap
from mylog import set_file_logger

# * configs
eps: float = 0.0001
min_r: float = 0.1
max_r: float = 2
step_r: float = 0.1
times_r: int = 10
nsample_sil: int = 50000
times_sil: int = 6

# * set debug or not
debug: bool = True if snakemake.params['debug'] > 0 else False

# FIXME: add cll and current i to filter barcodes
if not debug:
    logfnm: str = snakemake.log[0]
    b2g_fnm: str = snakemake.input['b2g_fnm']
    allen_k8: str = snakemake.input['allen_k8']
    ann_h5ad: str = snakemake.input['ann_h5ad']
    out_ann: str = snakemake.output['ann']
    out_cellmeta: str = snakemake.output['cellmeta']
    out_silht: str = snakemake.output['silht']
    npc: int = snakemake.params['npc']
    k_knn: int = snakemake.params['k_knn']
    knn_method: str = snakemake.params['knn_method']
    leiden_weight: float = snakemake.params['leiden_weight']
    threads: int = snakemake.threads
else:
    work_dir = f"{proj_dir}/01.clustering"
    out_dir = f'{work_dir}/out/test'
    os.makedirs(out_dir, exist_ok=True)
    logfnm = f'{out_dir}/test.clustering.log'
    b2g_fnm = f"{work_dir}/src/test/resource/b2g_test.csv"
    ann_h5ad = os.path.join(work_dir,
                            "src/test/resource",
                            "scRNAseq_sa2_test.ann.h5ad")
    allen_k8 = f"{proj_dir}/meta/AIT21_k8_markers.txt"
    npc = 50
    k_knn = 50
    knn_method = "kdtree"
    leiden_weight = 1
    threads = 2
    times_sil = 2
    times_r = 3
    min_r = 0.5
    max_r = 0.8
    step_r = 0.1
    prefix = f"RNA_k8_npc{npc}_k{k_knn}_L1"
    out_ann = os.path.join(out_dir, f"{prefix}_0.h5ad")
    out_cellmeta = os.path.join(out_dir, f"{prefix}_0.csv")
    out_silht = os.path.join(out_dir, f"{prefix}_0_silht.csv")

if knn_method == "pynndescent":
    raise RuntimeError(
        f"SnapATAC2 with knn method {knn_method} has bug.")

# * set logger
logger: Logger = set_file_logger(logfnm, name="sa2_clustering")

logger.info(f"SnapATAC2 version: {sa2.__version__}.")
if debug:
    logger.info("Under DEBUG mode.")
else:
    logger.info("Under NORMAL mode.")

# * load input
b2g: pd.DataFrame = pd.read_csv(b2g_fnm, sep=",", header=0)
with open(allen_k8, 'r') as f:
    k8s = [g.strip() for g in f.readlines()]
ann: sa2.AnnData = sa2.read(filename=ann_h5ad, backed='r')

# * use k8 features
# FIXME support scipy variable features
# this variable feature selection can be run on both full scale
# and on k8 features
# Let's run first round clustering, then for the next round, we
# can work on this.

genes_ann: List[str] = ann.var_names
k8_ann: List[str] = [True if g in k8s else False for g in genes_ann]
logger.info(f"After intersect with k8_allen, {sum(k8_ann)} features left.")
ann_var: sa2.AnnData = ann.subset(var_indices=k8_ann, out=out_ann)
ann.close()

# * select features with variations
logger.info("Select features with variations.")
ann_in_mem = ann_var.to_memory()
max_exp = ann_in_mem.X.max(axis=0).toarray()[0, :]
min_exp = ann_in_mem.X.min(axis=0).toarray()[0, :]
delta_exp = max_exp - min_exp
selected = polars.Series(delta_exp > eps)
ann_var.var['selected'] = selected
selected = delta_exp > eps
logger.info(f"{sum(selected)} features are kept.")
del ann_in_mem

# * non-linear embed
logger.info("Run spectral embedding.")
sa2.tl.spectral(adata=ann_var, n_comps=npc,
                features='selected',
                weighted_by_sd=True,
                inplace=True)
# * umap
logger.info("Run UMAP with given a and b.")
myumap.umap(adata=ann_var,
            n_comps=2,
            use_rep='X_spectral',
            key_added="umap",
            inplace=True,
            # use no seed for parallism
            # random_state=0,
            random_state=None,
            n_jobs=threads,
            a=1.8956,
            b=0.8006,
            init='spectral',
            metric='euclidean')
# * knn
logger.info(f"Run KNN with k={k_knn} and method={knn_method}.")
sa2.pp.knn(adata=ann_var,
           n_neighbors=k_knn,
           use_dims=None,
           use_rep='X_spectral',
           method=knn_method,
           inplace=True,
           random_state=0)

# * leiden and silhouette


def run_leiden(r: float) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Just a wrapper, it reads variables in the outside.
    """
    logger.info(f"Run Leiden under reso {r} with repeat {times_r}")
    leidens = np.ndarray(shape=(ann_var.shape[0], times_r), dtype='int')
    sils = np.ndarray(shape=(times_sil, times_r), dtype=float)
    for i in range(times_r):
        logger.info(f"Current randome state is {i}.")
        # dtype of element in leidens: '<U2'.
        t = sa2.tl.leiden(adata=ann_var,
                          resolution=r,
                          objective_function='modularity',
                          min_cluster_size=5,
                          random_state=i,
                          key_added='leiden',
                          use_leidenalg=False,
                          weighted=leiden_weight,
                          inplace=False)
        leidens[:, i] = t.astype(int)
        if np.unique(leidens[:, i]).shape[0] < 2:
            logger.info("No more clusters in current clustering.")
            logger.info("Then set silhouette as 0.0.")
            sils[:, i] = 0.0
            continue
        if ann_var.shape[0] < nsample_sil:
            logger.info("Get silhouette scores with all samples.")
            sils[:, i] = silhouette_score(X=ann_var.obsm['X_spectral'],
                                          labels=leidens[:, i],
                                          sample_size=None,
                                          random_state=None)
            continue
        logger.info(f"Get silhouette scores with repeat {times_sil}.")
        for j in range(times_sil):
            sils[j, i] = silhouette_score(X=ann_var.obsm['X_spectral'],
                                          labels=leidens[:, i],
                                          sample_size=nsample_sil,
                                          random_state=j)
    r1 = pd.DataFrame(data=leidens,
                      index=ann_var.obs_names,
                      columns=[f"r{r:.1f}_seed{j}" for j in range(times_r)])
    r2 = pd.DataFrame(data=sils,
                      index=[f"seed{j}" for j in range(times_sil)],
                      columns=r1.columns)
    return (r1, r2)


logger.info("Run leiden with resolutions.")
logger.info(f"From {min_r} to {max_r} with step {step_r}.")
rs: np.ndarray = np.arange(min_r, max_r+step_r, step_r)
# make sure all resos:0.1,0.2 .., 1, ..,2.
# Otherwise, some float might be 0.20000004
rs = np.array([round(i, 1) for i in rs])
leisil: List[Tuple[pd.DataFrame, pd.DataFrame]] = list(map(run_leiden, rs))
# logger.info(f"Use {threads} threads in parallel.")
# - not work in interactive mode
# - will generate ann_var again, which cause problem.
# with Pool(threads) as p:
#     leisil: List[Tuple[np.ndarray, np.ndarray]] = p.apply_async(
#         run_leiden, rs).get(9999999)
logger.info("Finish to save ann data with embedding, knn and umap.")
# * summarize results
barcodes = pd.DataFrame(data=ann_var.obs_names,
                        dtype=str,
                        index=ann_var.obs_names,
                        columns=['barcode'])
leidens: pd.DataFrame = pd.concat(objs=[i[0] for i in leisil],
                                  axis=1)
umap = pd.DataFrame(data=ann_var.obsm['X_umap'],
                    index=ann_var.obs_names,
                    columns=["umap1", "umap2"],
                    dtype=float)
cell_meta: pd.DataFrame = pd.concat([barcodes, leidens, umap], axis=1)
cell_meta.to_csv(out_cellmeta,
                 sep=',',
                 header=True, index=False)
logger.info("Finish to save cell_meta.")
# * calculate avg sils
sils: pd.DataFrame = pd.concat(
    [i[1] for i in leisil], axis=1).transpose()
sils['reso_seed'] = sils.index
sils.to_csv(out_silht, sep=',', header=True, index=False)
logger.info("Finish to save silhouette scores.")
ann_var.close()
