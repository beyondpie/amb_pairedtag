"""
ATAC-seq data are pre-handled under 05.CRE/src/main.python/sup.getgmat.sa2v26.snATAC.py
"""
import os
import pandas as pd
import anndata as ad

import tmpPypkg.globalvar as gv

# * meta
k8g = pd.read_csv(
    os.path.join(gv.pt_projd, "meta", "AIT21_k8_markers.txt"),
    header=None).iloc[:, 0]
snATACd = os.path.join(gv.pt_projd, "data/snATAC")

# * prepare nn cells
ann_nn = ad.read_h5ad(
    os.path.join(snATACd, "snATAC_gmat_obsv9.8_ptregion_nn_ann.h5ad"),
    backed=None)

k8g_ann = ann_nn.var_names.isin(k8g)
use_cells = ann_nn.obs.groupby(
    'subclass_id_label_v3', observed=True).apply(
    lambda x: x.sample(n=min(10000, x.shape[0]), random_state=0)
).index.droplevel(0)

ann_nnds = ann_nn[use_cells, ann_nn.var_names.isin(k8g)].copy()
ann_nnds.write_h5ad(
    os.path.join(snATACd, "snATAC.gmat.nn.ds.h5ad"))

# * prepare neu cells
ann_neu = ad.read_h5ad(
    os.path.join(snATACd, "snATAC_gmat_obsv9.8_ptregion_neu_ann.h5ad"),
    backed=None
)
k8g_ann = ann_neu.var_names.isin(k8g)
use_cells = ann_neu.obs.groupby(
    "subclass_id_label_v3", observed=True).apply(
        lambda x: x.sample(n=min(1000, x.shape[0]), random_state=0)
).index.droplevel(0)

ann_neuds = ann_neu[use_cells, ann_neu.var_names.isin(k8g)].copy()
ann_neuds.write_h5ad(
    os.path.join(snATACd, "snATAC.gmat.neu.ds.h5ad")
)

# * 2024-11-13
# prepare non-down-sample snATAC data
ann_nn = ad.read_h5ad(
    os.path.join(snATACd, "snATAC_gmat_obsv9.8_ptregion_nn_ann.h5ad"),
    backed=None)
ann_nn = ann_nn[:, ann_nn.var_names.isin(k8g)].copy()
ann_nn.write_h5ad(
    os.path.join(snATACd, "snATAC.ptregion.gmat.nn.h5ad")
)

ann_neu = ad.read_h5ad(
    os.path.join(snATACd, "snATAC_gmat_obsv9.8_ptregion_neu_ann.h5ad"),
    backed=None)
ann_neu = ann_neu[:, ann_neu.var_names.isin(k8g)].copy()
ann_neu.write_h5ad(
    os.path.join(snATACd, "snATAC.ptregion.gmat.neu.h5ad")
)

# * 2024-11-18
# merge downsampled nn and neu
ann_nnds = ad.read_h5ad(os.path.join(snATACd, "snATAC.gmat.nn.ds.h5ad"),
                        backed=None)
ann_neuds = ad.read_h5ad(os.path.join(snATACd, "snATAC.gmat.neu.ds.h5ad"),
                         backed=None)
ann_ds = ad.concat([ann_nnds, ann_neuds], merge="same")
ann_ds.write_h5ad(os.path.join(snATACd, "snATAC.ptregion.gmat.k8.ds.raw.h5ad"))

