import argparse
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from pathlib import Path
from harmonypy import run_harmony
# from snapatac2 python-contrib path
import integration as it
# from local python package
from tmpPypkg.utils import set_file_logger
from tmpPypkg.singlecell import normalize_data, scanpy_PCA_plus
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(
  prog="run.cca",
  description="run python CCA."
)
parser.add_argument(
  '--ref', type=str, help="scanpy AnnData as ref")
parser.add_argument(
  '--refnm', type=str, default="ref", help="batch name for ref")
parser.add_argument(
  '--query', type=str, help="scanpy AnnData as query")
parser.add_argument(
  '--querynm', type=str, default="query", help="batch name for query")


# normalization
parser.add_argument("--refnorm", action="store_true",
                    help="perform normalization on ref.")
parser.add_argument("--querynorm", action="store_true",
                    help="perform normalization on query.")


# feature selection
parser.add_argument(
    "--reflavor", type=str, default="seurat_v3",
    help="feature selection flavor for reference")
parser.add_argument(
    "--queryflavor", type=str, default="seurat_v3",
    help="feature selection flavor for query")
parser.add_argument(
  '--refhv', action='store_true', help="use ref highly variable features")
parser.add_argument(
  '--queryhv', action='store_true', help="use query highly variable features")
parser.add_argument(
  '--nhv', type=int, default=3000, help="# highly variable features")

# PCA
parser.add_argument(
    "--refPCA", action="store_true", help="run PCA for ref")
parser.add_argument(
    "--queryPCA", action="store_true", help="run PCA for query")
parser.add_argument(
    "--npca", type=int, default=50, help="#components for PCA")

# integration with CCA
parser.add_argument(
  '--k_anchor', type=int, default=5, help="#neighbors for picking anchors")
parser.add_argument(
  '--k_filter', type=int, default=200, help="#neighbors for filtering anchors")
parser.add_argument(
  '--k_score', type=int, default=30, help="#neighbors for scoring anchors")
parser.add_argument(
  '--ncca', type=int, default=50, help="#components for CCA")
parser.add_argument(
  '--maxfea', type=int, default=200, help="max #features for filtering anchors")
parser.add_argument(
  '--maxcell', type=int, default=200000, help="max #cells used in exact cca")
parser.add_argument(
  '--k_local', type=int, default=0, help="#neighbors to adjust anchor scores")
parser.add_argument(
  '--k_weight', type=int, default=100, help="#neighbors for weighting anchors")
parser.add_argument(
  '--anchorpath', type=str, help="direcotry for saving anchors")
parser.add_argument(
  '--inth5ad', type=str, help="scanpy AnnData fnm for saving results")
parser.add_argument(
  '--logfnm', type=str, help="logfnm for recording process")
parser.add_argument(
  '--flagfnm', type=str, help="flagfnm to record job finish status")
parser.add_argument(
  '--run_harmony', action='store_true', help="run harmony")
parser.add_argument(
  '--run_umap', action="store_true", help="run UMAP")
parser.add_argument(
    '--umapfnm', type=str,
    default="test.pdf",
    help="figure name to save umap")

# * parse args
args = parser.parse_args()
# DEBUG only
# use list inside parse_args for debugs
# print("DEBUG...")
# projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
# datad = f"{projd}/data"
# outd = f"{projd}/03.integration/out"
# args = parser.parse_args(
#     ["--ref", f"{datad}/allen/allen.10xv3.pt.nn.ds.h5ad",
#      "--refnm", "allen_nn",
#      "--query", f"{datad}/snmC_snm3C/snmC.scaled.gmat.nn.ds.h5ad",
#      "--querynm", "snmC_nn",
#      "--refnorm",
#      "--reflavor", "seurat_v3", "--refhv",
#      "--queryflavor", "seurat", "--queryhv",
#      "--refPCA", "--queryPCA",
#      "--anchorpath", f"{outd}/int_allen_snmC_nn",
#      "--inth5ad", f"{outd}/int_allen_snmC_nn/intgn.h5ad",
#      "--logfnm", f"{outd}/int_allen_snmC_nn/log.txt",
#      "--flagfnm", f"{outd}/int_allen_snmC_nn/flag.done",
#      "--run_harmony", "--run_umap",
#      "--umapfnm", f"{outd}/int_allen_snmC_nn/umap.pdf"
#      ]
# )


log = set_file_logger(
  fnm=args.logfnm, fmode='a', name='Run CCA-based integration',
  log_level=logging.DEBUG)

# * load data
ref = ad.read_h5ad(args.ref, backed=None)
query = ad.read_h5ad(args.query, backed=None)

# * normalization
# CPM -> log1p -> scale
if args.refnorm:
    log.info("Perform normalization for ref.")
    normalize_data(ref)

if args.querynorm:
    log.info("Perform normalization for query.")
    normalize_data(query)

# * add single-cell feature selection
# here obs_names is pd.Index
init_genes = ref.var_names.intersection(query.var_names)
reflavor = args.reflavor
queryflavor = args.queryflavor
log.info(
  f"Use all the {init_genes.shape[0]} joint feature as init features.")

reflayer = "raw_counts" if args.reflavor == "seurat_v3" else None
reflayer = "raw_counts" if "raw_counts" in ref.layers else None

querylayer = "raw_counts" if args.queryflavor == "seurat_v3" else None
querylayer = "raw_counts" if "raw_counts" in query.layers else None

if args.refhv and (not args.queryhv):
    df = sc.pp.highly_variable_genes(
      adata=ref, layer=reflayer, n_top_genes=args.nhv, flavor=reflavor,
      inplace=False)
    init_genes = ref.var_names[df.highly_variable]
    log.info(
      f"Use ref {init_genes.shape[0]} variable genes.")
elif (not args.refhv) and (args.queryhv):
    df = sc.pp.highly_variable_genes(
      adata=query, layer=querylayer, n_top_genes=args.nhv, flavor=queryflavor,
      inplace=False)
    init_genes = query.var_names[df.highly_variable]
    log.info(
      f"Use query {init_genes.shape[0]} variable genes.")
elif (args.refhv) and (args.queryhv):
    df1 = sc.pp.highly_variable_genes(
      adata=ref, layer=reflayer, n_top_genes=args.nhv, flavor=reflavor,
      inplace=False)
    df2 = sc.pp.highly_variable_genes(
      adata=query, layer=querylayer, n_top_genes=args.nhv, flavor=queryflavor,
      inplace=False)
    gene_names_ref = ref.var_names[df1.highly_variable]
    gene_names_query = query.var_names[df2.highly_variable]
    init_genes = gene_names_ref.intersection(gene_names_query)
    log.info(
      f"Use jointly {init_genes.shape[0]} variable genes.")
else:
    log.info(f"Use all the {init_genes.shape[0]} genes.")

ref = ref[:, init_genes].copy()
query = query[:, init_genes].copy()

# * run PCA
if args.refPCA:
    log.info("Run PCA for reference.")
    scanpy_PCA_plus(ref, n_comps=args.npca,
                    weight_by_var=True, zero_center=True)
if args.queryPCA:
    log.info("Run PCA for query.")
    scanpy_PCA_plus(query, n_comps=args.npca,
                    weight_by_var=True, zero_center=True)


# * find anchor
log.info("Start to find anchor.")
si = it.SeuratIntegration(n_jobs=4, random_state=0)
si.find_anchor(
    adata_list=[ref, query],
    k_local=None if args.k_local < 1 else args.k_local,
    k_anchor=args.k_anchor,
    key_anchor="X",
    dim_red="cca",
    max_cc_cells=args.maxcell,
    k_score=args.k_score,
    k_filter=args.k_filter,
    scale1=True,
    scale2=True,
    n_components=args.ncca,
    n_features=args.maxfea
)


# * run integration
log.info("Perform integration.")
intgrn = si.integrate(
  key_correct="X_pca",
  row_normalize=True,
  k_weight=args.k_weight,
  sd=1,
  alignments=[[[0], [1]]]
)

ref_obs = pd.DataFrame(
  {
    "barcode": ref.obs_names,
    "batch": args.refnm
  },
  index=ref.obs_names
)
query_obs = pd.DataFrame(
  {
    "barcode": query.obs_names,
    "batch": args.querynm
  },
  index=query.obs_names
)

# * organize the results
ann_intgrn = ad.AnnData(
  X=np.concatenate(intgrn, axis=0),
  obs=pd.concat([ref_obs, query_obs]),
  obsm={
    'raw_pca':
    np.concatenate([ref.obsm['X_pca'], query.obsm['X_pca']], axis=0)
  }
)

# * perform harmony
if args.run_harmony:
    log.info("Perform Harmony.")
    hm = run_harmony(ann_intgrn.X, ann_intgrn.obs, 'batch',
                     max_iter_harmony=30).Z_corr.T
    log.info("Finish Harmony.")
    ann_intgrn.obsm['harmony'] = hm

# * perrorm UMAP
if args.run_umap:
    log.info("Perfom UMAP.")
    sc.pp.neighbors(
      ann_intgrn, method="umap", metric="euclidean",
      use_rep="X")
    sc.tl.leiden(ann_intgrn, resolution=0.5)
    try:
        log.info("Try init the param using paga.")
        sc.tl.paga(ann_intgrn, groups='leiden')
        sc.pl.paga(ann_intgrn, plot=False)
        sc.tl.umap(ann_intgrn, min_dist=0.1, init_pos='paga')
    except Exception:
        log.info(
          'Init with PAGA failed, use default spectral init')
        sc.tl.umap(ann_intgrn, min_dist=0.1)
    finally:
        log.info("Finish UMAP.")

ann_intgrn.write_h5ad(args.inth5ad)
log.info(f"Saved result to {args.inth5ad}.")
log.info(f"Saved anchor object to {args.anchorpath}.")

# Error when setting all the members as True
# TypeError: Cannot setitem on a Categorical with a new category (nan), set the categories first
si.save(args.anchorpath,
        save_local_knn=False,
        save_raw_anchor=False,
        save_mutual_knn=False,
        save_adata=False)
Path(args.flagfnm).touch()

if args.run_umap:
    log.info(f"Generate UMAP figure to {args.umapfnm}.")
    fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
    d = pd.DataFrame(
        {
            "x": ann_intgrn.obsm['X_umap'][:, 0],
            "y": ann_intgrn.obsm['X_umap'][:, 1],
            "batch": ann_intgrn.obs['batch']
        },
        index=ann_intgrn.obs.index
    )

    sns.scatterplot(
        x="x",
        y="y",
        data=d,
        ax=ax,
        hue="batch"
    )
    fig.savefig(args.umapfnm)

log.info("Jobs done, and good luck!")
