import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import snapatac2 as sa2
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
from typing import Dict


# from snapatac2 python-contrib path
import integration as it
from integration.seurat_class import scanpy_PCA_plus
# from local python package
from tmpPypkg import globalvar as gv
from tmpPypkg import utils
# utils.reimport("integration")
from harmonypy import run_harmony

# * functions
def normalize_data_and_pca(adata: ad.AnnData) -> None:
    print("Save raw counts into layer.")
    adata.layers['raw_counts'] = adata.X
    print("Perform logCPM normalization for raw_counts.")
    sc.pp.normalize_total(adata, target_sum=1e6,
                          exclude_highly_expressed=False,
                          max_fraction=0.05,
                          key_added=None,
                          inplace=True,
                          layer="raw_counts",
                          copy=False)
    sc.pp.log1p(adata, base=np.e)
    print("Scaling the logCPM.")
    sc.pp.scale(adata, zero_center=True)
    print("Perform PCA.")
    scanpy_PCA_plus(
        adata, n_comps=50, weight_by_var=True,
        zero_center=True)
    return None


# * prepare snapatac2 object for integration
# -  nn
allenNNscann_fnm = os.path.join(
    gv.pt_projd,
    "data/",
    ".10xv3.pt.regions.nn.imn.ann.k8.rawcount.h5ad"
)
a = ad.read_h5ad(allenNNscann_fnm, backed=None)

# - pt nn
ptNNsa2ann_fnm = os.path.join(
    gv.pt_projd,
    "data/pairedtag_ann",
    "pt.RNAseq.sa2ann.NN.k8.h5ad"
)
# once loading into memory, it's AnnData even the obs field.
p = sa2.read(filename=ptNNsa2ann_fnm, backed=None)
# add cellmeta to pairedtag data
ptcm = pd.read_csv(gv.ptcellmetafnm,
                   sep=",", header=0)
ptcm.set_index(keys='barcode', drop=False, inplace=True)
old_obs = p.obs
new_obs = ptcm.loc[p.obs.index]
p.obs = new_obs

# perform downsample and save them
# for paired-tag data
pnds = 30
p_obs_ds = new_obs.groupby(
  'pairedtagCluster').apply(
    lambda x: x.sample(n=min(pnds, x.shape[0]), random_state=0))
pp = p[p_obs_ds.barcode, ].copy()
# remove columns with nan values
old_obs = pp.obs
new_obs = old_obs.copy()
new_obs.drop(columns=['rannot.cl', 'rannot.sp'], inplace=True)
pp.obs = new_obs
# write to file
pp.write_h5ad(
  filename=os.path.join(gv.pt_projd, "data/pairedtag_ann",
                        "pt.RNAseq.ann.NN.k8.nds30.h5ad"))
# for  dataset
a_obs = a.obs.copy()
a_obs.insert(0, 'barcode', a_obs.index)
ands = 1000
a_obs_ds = a_obs.groupby(
  "cl").apply(
    lambda x: x.sample(n=min(ands, x.shape[0]), random_state=0))
aa = a[a_obs_ds.barcode, ].copy()
aa.write_h5ad(
  filename=os.path.join(gv.pt_projd, "data/",
    ".10xv3.pt.regions.nn.imn.ann.k8.rawcount.nds1000.h5ad")
)

# 2. CCA
a = ad.read_h5ad(
    os.path.join(
        gv.pt_projd,
        "data/",
        ".10xv3.pt.regions.nn.imn.ann.k8.rawcount.nds1000.h5ad")
)

p = ad.read_h5ad(
    os.path.join(
        gv.pt_projd, "data/pairedtag_ann",
        "pt.RNAseq.ann.NN.k8.nds30.h5ad")
)

# perform PCA with scaling
gene_names = [g for g in p.var_names if g in a.var_names]
aa = a[:, gene_names].copy()
# aa.obsm['batch'] = np.array(["Allen-10Xv3-scRNA"] * aa.shape[0])
aa.obs['batch'] = 'Allen-10Xv3-scRNA'
normalize_data_and_pca(adata=aa)

pp = p[:, gene_names].copy()
pp.obs['batch'] = 'PairedTag-snRNA'
# pp.obsm['batch'] = np.array(["PairedTag-snRNA"] * pp.shape[0])
normalize_data_and_pca(adata=pp)

aa.write("test_integration/adata/aa.h5ad")
pp.write("test_integration/adata/pp.h5ad")

si = it.SeuratIntegration(n_jobs=4, random_state=0)
si.find_anchor(
    adata_list=[aa, pp],
    k_local=None,
    k_anchor=5,
    key_anchor="X",
    dim_red="cca",
    max_cc_cells=100000,
    k_score=30,
    k_filter=200,
    scale1=True,
    scale2=True,
    n_components=50,
    n_features=200
)

# when saving with keeping the attributes, nan cannot be saved error.
si.save("test_integration")
# si.load("test_integration")

corrected = si.integrate(
    key_correct="X_pca",
    row_normalize=True,
    k_weight=100,
    sd=1,
    alignments=[[[0], [1]]],
)

# adata_merge.obsm["X_pca_integrate"] = np.concatenate(corrected)
x_pca_int = np.concatenate(corrected)
np.save(file="test_integration/x_pca_int.npy", arr=x_pca_int)

# run harmony on x_pca_int
aa = ad.read_h5ad("test_integration/adata/aa.h5ad")
pp = ad.read_h5ad("test_integration/adata/pp.h5ad")
x_pca_int = np.load("test_integration/x_pca_int.npy")
aa.obs['barcode'] = aa.obs.index
aa.obs['batch'] = "Allen-10x-RNA"
pp.obs['batch'] = "pairedtag-RNA"
obs_meta = pd.concat([aa.obs[['barcode', 'batch']],
                      pp.obs[['barcode', 'batch']]])
obs_meta.set_index('barcode', drop=False, inplace=True)
pca_int = run_harmony(x_pca_int, obs_meta, 'batch',
                      max_iter_harmony=30).Z_corr.T

adata_merge = aa.concatenate(pp, batch_categories=["", "paireddtag"], batch_key="merged_batch",
                             index_unique=None)
aa.obs = aa.obs[['barcode', 'batch']]
pp.obs = pp.obs[['barcode', 'batch']]
a3 = aa.copy()
p3 = pp.copy()

adata_merge = ad.AnnData(
  X=np.concatenate([aa.obsm['X_pca'], pp.obsm['X_pca']], axis=0),
  obs=obs_meta,
  obsm={"pca_int": pca_int}
)
# sc.pp.neighbors(adata_merge, method="umap", use_rep="pca_int")
sc.tl.umap(adata_merge, min_dist=0.5,
           init_pos='spectral',
           alpha=1.0,
           gamma=1.0,
           a=None,
           b=None)
adata_merge.write_h5ad("test_integration/adata/adata_merge.h5ad")

fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
d = pd.DataFrame(
  {
    "x": adata_merge.obsm['X_umap'][:, 0],
    "y": adata_merge.obsm['X_umap'][:, 1],
    "batch": adata_merge.obs['batch']
  },
  index=adata_merge.obs.index
)

sns.scatterplot(
  x="x",
  y="y",
  data=d,
  ax=ax,
  hue="batch"
)
fig.savefig("test_NN_int.png")

# * now use this to co-emebding ATAC data
aa = ad.read_h5ad("test_integration/adata/aa.h5ad")

# load snATAC data and perform downsample
sn = ad.read_h5ad(
    os.path.join(gv.pt_projd, "data/snATAC",
                 "snATAC_gmat_obsv9.8_ptregion_nn_ann.h5ad"),
    backed=None)
gene_names = aa.var_names.intersection(sn.var_names)

sn_obs = sn.obs.copy()
sn_obs.insert(0, 'ubarcode', sn_obs.index)
sn_obs_ds = sn_obs.groupby('subclass_id_label_v3').apply(
    lambda x: x.sample(n=min(2000, x.shape[0]), random_state=0)
).copy()

snds = sn[sn_obs_ds.ubarcode, gene_names].copy()
snds.obs['batch'] = "snATAC"
normalize_data_and_pca(adata=snds)
snds.write("test_integration/adata/snds_nn.h5ad")
aa2 = aa[:, gene_names].copy()

si2 = it.SeuratIntegration(n_jobs=4, random_state=0)
# will use all the cpu during CCA,
# which is out of control of n_jobs
si2.find_anchor(
    adata_list=[aa2, snds],
    k_local=None,
    k_anchor=5,
    key_anchor="X",
    dim_red="cca",
    max_cc_cells=100000,
    k_score=30,
    k_filter=200,
    scale1=True,
    scale2=True,
    n_components=50,
    n_features=200
)

corct2 = si2.integrate(
    key_correct="X_pca",
    row_normalize=True,
    k_weight=100,
    sd=1,
    alignments=[[[0], [1]]],
)

corct_allen_pt = np.load("test_integration/x_pca_int.npy")
corct_al_pt_at = np.concatenate([corct_allen_pt, corct2[1]])
pp = ad.read_h5ad("test_integration/adata/pp.h5ad")

aa_obs = aa.obs.copy()
pp_obs = pp.obs.copy()
at_obs = snds.obs.copy()

aa_obs['barcode'] = aa.obs.index
aa_obs['batch'] = "allen"
pp_obs['batch'] = "pt"
at_obs['barcode'] = at_obs.index
at_obs['batch'] = 'atac'


obs_meta = pd.concat([
    aa_obs[['barcode', 'batch']],
    pp_obs[['barcode', 'batch']],
    at_obs[['barcode', 'batch']]
])

hm = run_harmony(corct_al_pt_at, obs_meta, 'batch',
                 max_iter_harmony=30).Z_corr.T

adata_merge = ad.AnnData(
  X=np.concatenate(
      [aa.obsm['X_pca'], pp.obsm['X_pca'], snds.obsm['X_pca']], axis=0),
  obs=obs_meta,
  obsm={"pca_int": hm}
)
sc.pp.neighbors(adata_merge, use_rep="pca_int")
resolution = 0.5
sc.tl.leiden(adata_merge, resolution=resolution)
# 2 hours to run 4M + 3K cell
# 4 hours to run 4M + 3K cell if using spectral init, the init step is very slow
min_dist = max(0.1, 1 - adata_merge.shape[0] / 60000)
try:
    sc.tl.paga(adata_merge, groups='leiden')
    sc.pl.paga(adata_merge, plot=False)
    sc.tl.umap(adata_merge, min_dist=min_dist, init_pos='paga')
except Exception:
    print('Init with PAGA failed, use default spectral init')
    sc.tl.umap(adata_merge, min_dist=min_dist)
adata_merge.write_h5ad(
    "test_integration/adata/adata_merge_al_pt_at.h5ad")


fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
d = pd.DataFrame(
  {
    "x": adata_merge.obsm['X_umap'][:, 0],
    "y": adata_merge.obsm['X_umap'][:, 1],
    "batch": adata_merge.obs['batch']
  },
  index=adata_merge.obs.index
)

sns.scatterplot(
  x="x",
  y="y",
  data=d,
  ax=ax,
  hue="batch"
)
fig.savefig("test_NN_int_al_pt_at.png")

adata = ad.read_h5ad("test_integration/adata/adata_merge_al_pt_at.h5ad")

# * test reload si and perform intergration
def simply_load_SeuratIntegration(
    model_path: str,
    adata_dict: Dict[int, ad.AnnData]
) -> it.SeuratIntegration:
    si = joblib.load(model_path)
    si.adata_dict = adata_dict
    return si


outd = os.path.join(gv.pt_projd, "03.integration",
                    "out/int_allen_snmC_nn")
adata_intgn = ad.read_h5ad(f"{outd}/intgn.h5ad", backed=None)
adata_ref = adata_intgn[adata_intgn.obs.batch == "allen_nn", ].copy()
adata_query = adata_intgn[adata_intgn.obs.batch == "snmC_nn", ].copy()

si = simply_load_SeuratIntegration(
    model_path=f"{outd}/model.lib",
    adata_dict={
        0: adata_ref,
        1: adata_query
    }
)

redo_intgn = si.integrate(
    key_correct="raw_pca",
    row_normalize=True,
    k_weight=100,
    sd=1,
    alignments=[[[0], [1]]]
)
