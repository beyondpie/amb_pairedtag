import sys
import os
import re
from typing import Dict
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from harmonypy import run_harmony

import integration as it
import tmpPypkg.globalvar as gv
from tmpPypkg.singlecell import normalize_data


# * meta
anchord = os.path.join(gv.pt_projd, "03.integration", "out")
annd = os.path.join(gv.pt_projd, "data")
ntop = 7500
npca = 50
# prefix = "int"
# suffix = "ds"
# coembed_class = "neu"
# allen_ann_fnm = os.path.join(
#     annd, "allen",
#     f"allen.10xv3.pt.{coembed_class}.ds.h5ad")
# pt_ann_fnm = os.path.join(
#     annd, "pairedtag_ann",
#     f"pt.RNA.ann.{coembed_class}.neu.withimn.k8.ds.h5ad")
# snmC_ann_fnm = os.path.join(
#     annd, "snmC_snm3C", f"snmC.scaled.gmat.{coembed_class}.ds.h5ad"
# )
# atac_ann_fnm = os.path.join(
#     annd, "snATAC", f"snATAC.gmat.{coembed_class}.ds.h5ad")

# prefix = "all"
# suffix = "all"
# coembed_class = "neu"
# allen_ann_fnm = os.path.join(
#     annd, "allen",
#     f"allen.10xv3.pt.regions.{coembed_class}.imn.ann.k8.rawcount.h5ad")
# pt_ann_fnm = os.path.join(
#     annd, "pairedtag_ann", f"pt.RNA.ann.{coembed_class}.withimn.k8.h5ad")
# snmC_ann_fnm = os.path.join(
#     annd, "snmC_snm3C", f"snmC.scaled.gmat.k8.{coembed_class}.h5ad")
# atac_ann_fnm = os.path.join(
#     annd, "snATAC", f"snATAC.ptregion.gmat.{coembed_class}.h5ad")

# model path of seurat integrations
# si_pt_modelpath = os.path.join(anchord,
#                                 f"{prefix}_allen_pt_{coembed_class}", "model.lib")
# si_snmC_modelpath = os.path.join(anchord,
#     f"{prefix}_allen_snmC_{coembed_class}", "model.lib")
# si_atac_modelpath = os.path.join(anchord,
#                                  f"{prefix}_allen_snATAC_{coembed_class}",
#                                  "model.lib")

# prefix = "nneuds"
# prefix = "nneudspt"
prefix = sys.argv[1]
print(f"coembed for {prefix}.")
allen_ann_fnm = os.path.join(annd, "allen", "allen.10xv3.pt.k8.ds.raw.h5ad")
pt_ann_fnm = os.path.join(annd, "pairedtag_ann", "pt.RNA.ann.k8.ds.raw.h5ad")
snmC_ann_fnm = os.path.join(annd, "snmC_snm3C", "snmC.scaled.gmat.k8.ds.h5ad")
atac_ann_fnm = os.path.join(annd, "snATAC", "snATAC.ptregion.gmat.k8.ds.raw.h5ad")

si_pt_modelpath = os.path.join(anchord, f"{prefix}_allen_pt", "model.lib")
si_snmC_modelpath = os.path.join(anchord, f"{prefix}_allen_snmC", "model.lib")
si_atac_modelpath = os.path.join(anchord, f"{prefix}_allen_snATAC", "model.lib")


# * functions
def loadsi(model_path: str, adata_dict: Dict[int, ad.AnnData]) -> it.SeuratIntegration:
    import joblib

    si = joblib.load(model_path)
    si.adata_dict = adata_dict
    return si


def runUMAP(ann: ad.AnnData, use_rep: str = "X") -> ad.AnnData:
    print(f"Run UMAP with {use_rep} rep.")
    sc.pp.neighbors(ann_intgrn, method="umap", metric="euclidean", use_rep=use_rep)
    sc.tl.leiden(ann_intgrn, resolution=0.5)
    try:
        sc.tl.paga(ann_intgrn, groups="leiden")
        sc.pl.paga(ann_intgrn, plot=False)
        r = sc.tl.umap(ann_intgrn, min_dist=0.1, init_pos="paga", copy=True)
    except Exception:
        print("Init with PAGA failed, use default spectral init")
        r = sc.tl.umap(ann_intgrn, min_dist=0.1, copy=True)
    finally:
        print("Finish UMAP.")
    return r


# ** main **
# * load all the anndatas
print("Load all the anndatas.")
ann_allen = ad.read_h5ad(allen_ann_fnm, backed=None)
ann_pt = ad.read_h5ad(pt_ann_fnm, backed=None)
ann_snmC = ad.read_h5ad(snmC_ann_fnm, backed=None)
ann_atac = ad.read_h5ad(atac_ann_fnm, backed=None)

# * load all the Seurat Integrations
print("Load all the seurat integration results.")
si_pt = loadsi(model_path=si_pt_modelpath, adata_dict={0: ann_pt, 1: ann_allen})
si_snmC = loadsi(model_path=si_snmC_modelpath, adata_dict={0: ann_pt, 1: ann_snmC})
si_atac = loadsi(model_path=si_atac_modelpath, adata_dict={0: ann_pt, 1: ann_atac})

# * for each anndata, run PCA
# top variable features then set the joint features
print("Normalize data one by one.")
normalize_data(ann_allen)
normalize_data(ann_pt)
normalize_data(ann_atac)

print("Select highly variable genes.")
sc.pp.highly_variable_genes(
    ann_allen, layer="raw_counts", n_top_genes=ntop, flavor="seurat_v3"
)
sc.pp.highly_variable_genes(
    ann_pt, layer="raw_counts", n_top_genes=ntop, flavor="seurat_v3"
)
sc.pp.highly_variable_genes(ann_snmC, n_top_genes=ntop, flavor="seurat")
sc.pp.highly_variable_genes(
    ann_atac, layer="raw_counts", n_top_genes=ntop, flavor="seurat_v3"
)

allen_hv = ann_allen.var.index[ann_allen.var.highly_variable]
pt_hv = ann_pt.var.index[ann_pt.var.highly_variable]
snmC_hv = ann_snmC.var.index[ann_snmC.var.highly_variable]
atac_hv = ann_atac.var.index[ann_atac.var.highly_variable]

hvg = allen_hv[allen_hv.isin(pt_hv)]
hvg = hvg[hvg.isin(snmC_hv)]
# 1597 variable genes for downsample dataset.
# 6166 genes for nneu downsample
hvg = hvg[hvg.isin(atac_hv)]
print(f"{len(hvg)} genes are used.")

# run PCA for each modality
print("Run PCA for each modality.")
sc.pp.pca(ann_allen, n_comps=50, mask_var=ann_allen.var_names.isin(hvg))
sc.pp.pca(ann_pt, n_comps=50, mask_var=ann_pt.var_names.isin(hvg))
sc.pp.pca(ann_snmC, n_comps=50, mask_var=ann_snmC.var_names.isin(hvg))
sc.pp.pca(ann_atac, n_comps=50, mask_var=ann_atac.var_names.isin(hvg))

# * centered on Allen's RNA-seq, perform co-embed
print("Perform integration centered on Allen dataset.")
int_allen_pt = si_pt.integrate(
    key_correct="X_pca", row_normalize=True, k_weight=100, sd=1, alignments=[[[0], [1]]]
)

int_allen_snmC = si_snmC.integrate(
    key_correct="X_pca", row_normalize=True, k_weight=100, sd=1, alignments=[[[0], [1]]]
)

int_allen_atac = si_atac.integrate(
    key_correct="X_pca", row_normalize=True, k_weight=100, sd=1, alignments=[[[0], [1]]]
)

# * organize the int results
pt_pca = int_allen_pt[0]
allen_pca = int_allen_pt[1]
snmC_pca = int_allen_snmC[1]
atac_pca = int_allen_atac[1]
X = np.concatenate([allen_pca, pt_pca, snmC_pca, atac_pca], axis=0)

allen_obs = pd.DataFrame(
    {"barcode": ann_allen.obs_names, "batch": "Allen10Xv3"}, index=ann_allen.obs_names
)

pt_obs = pd.DataFrame(
    {"barcode": ann_pt.obs_names, "batch": "PariedTag"}, index=ann_pt.obs_names
)

snmC_obs = pd.DataFrame(
    {"barcode": ann_snmC.obs_names, "batch": "snmC"}, index=ann_snmC.obs_names
)

atac_obs = pd.DataFrame(
    {"barcode": ann_atac.obs_names, "batch": "snATAC"}, index=ann_atac.obs_names
)
obs = pd.concat([allen_obs, pt_obs, snmC_obs, atac_obs])

ann_intgrn = ad.AnnData(
    X=X,
    obs=obs,
    obsm={
        "raw_pca": np.concatenate(
            [
                ann_allen.obsm["X_pca"],
                ann_pt.obsm["X_pca"],
                ann_snmC.obsm["X_pca"],
                ann_atac.obsm["X_pca"],
            ],
            axis=0,
        )
    },
)

# * run harmony
print("Run Harmony")
hm = run_harmony(ann_intgrn.X, ann_intgrn.obs, "batch",
                 max_iter_harmony=30).Z_corr.T
ann_intgrn.obsm["X_harmony"] = hm

# * UMAP without harmony
print("Run UMAP without harmony.")
ann1 = runUMAP(ann=ann_intgrn, use_rep="X")
ann_intgrn.obsm["X_umap"] = ann1.obsm["X_umap"]

# * UMAP with harmony
print("Run UMAP with harmony.")
ann2 = runUMAP(ann=ann_intgrn, use_rep="X_harmony")
ann_intgrn.obsm["X_umap_hm"] = ann2.obsm["X_umap"]

# add class info for ann_intgrn
allenClMeta = pd.read_csv(
    f"{gv.pt_projd}/meta/AIT21_annotation_freeze_081523.tsv", sep="\t", header=0
)
sc2class = allenClMeta[["subclass_id_label", "class_id_label"]].drop_duplicates()
subclass_label = sc2class["subclass_id_label"].apply(lambda x: re.sub(r"^\d+ ", "", x))
sc2class.insert(0, "subclass_label", subclass_label)
sc2class.set_index("subclass_id_label", drop=False, inplace=True)
sclabel2class = sc2class.set_index("subclass_label", drop=False, inplace=False)


class_id_label = pd.Series(sc2class.class_id_label.unique())
class_label = class_id_label.apply(lambda x: re.sub(r"^\d+ ", "", x))

allen_barcode2class = pd.DataFrame(
    {
        "barcode": ann_allen.obs_names,
        "class": sc2class.loc[ann_allen.obs.subclass_id_label, "class_id_label"],
    }
)
allen_barcode2class.set_index("barcode", drop=False, inplace=True)

pt_barcode2class = pd.DataFrame(
    {
        "barcode": ann_pt.obs_names,
        "class": sc2class.loc[ann_pt.obs.subclass, "class_id_label"],
    }
)
pt_barcode2class.set_index("barcode", drop=False, inplace=True)

snmC_barcode2class = pd.DataFrame(
    {
        "barcode": ann_snmC.obs_names,
        "class": sclabel2class.loc[ann_snmC.obs.SubClass, "class_id_label"],
    }
)
snmC_barcode2class.set_index("barcode", drop=False, inplace=True)

atac_barcode2class = pd.DataFrame(
    {"barcode": ann_atac.obs_names, "class": ann_atac.obs["class_id_label_v3"]}
)

barcode2class = pd.concat(
    [allen_barcode2class, pt_barcode2class, snmC_barcode2class, atac_barcode2class]
)

ann_intgrn.obs.insert(3, "class", barcode2class.loc[ann_intgrn.obs_names, "class"])

ann_intgrn.write_h5ad(
    os.path.join(anchord, "coembed", f"{prefix}.ann.intgrn.umap.class.h5ad")
)

# pt_class = pt_barcode2class["class"].unique()
# ann_intgrn = ann_intgrn[ann_intgrn.obs["class"].isin(pt_class),]
print("Done. Good luck.")
