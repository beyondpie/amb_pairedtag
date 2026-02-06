"""
Perform SnapATAC2 clustering for PairedTag dataset.
- currently for all the pairedtag data.
"""
import os
import numpy as np
import pandas as pd
import pyprojroot
import scanpy as sc
import anndata
import snapatac2 as sa2

# * configs
proj_root = pyprojroot.here()
meta_dir = os.path.join(proj_root, "meta")
out_dir = os.path.join(proj_root,
                       "01.clustering", "out")

# * parameters
npcs = 50
knn_k = 50
knn_method = "annoy"
reso = 0.8

# * load data
ann = sa2.read(
    filename = os.path.join(out_dir, "scRNAseq_sa2_all.ann.h5ad"),
    backed = None
)

# * load allen's features
allen_k8_genes = pd.read_csv(
    os.path.join(meta_dir, "AIT21_k8_markers.txt"),
    sep = ",", header = True)
allen_genes = allen_k8_genes[0].to_list()

# * subset ann
ann_k8 = ann[:, ann.var_names.isin(allen_genes)]

# * spectral embedding
# about 20 minutes
sa2.tl.spectral(
    adata = ann, n_comps = npcs,
    features = None,
    weighted_by_sd = True,
    inplace = True
)

# umap
sa2.tl.umap(
    adata = ann,
    n_comps = 2,
    use_rep = 'X_spectral',
    key_added = 'umap',
    inplace = True
)

# * knn
sa2.pp.knn(
    adata = ann,
    n_neighbors = knn_k,
    use_dims = None,
    use_rep = 'X_spectral',
    method = knn_method,
    inplace = True,
    random_state = 0
)

# * clustering
sa2.tl.leiden(
    adata = ann,
    resolution = reso,
    objective_function = 'modularity',
    random_state = 0,
    key_added = 'leiden',
    use_leidenalg = False,
    weighted = False,
    inplace = True
)

# visualize
outf_umap = os.path.join(
    out_dir, "figures", "scRNA_L1_UMAP_k8_50nc_umap.pdf")
umap_plt = sa2.pl.umap(adata = ann,
                       show = False,
                       color = 'leiden',
                       interactive = False,
                       out_file = outf_umap)
# * output
l1_out_dir = os.path.join(
    out_dir, "test_L1_sa2")
os.makedirs(l1_out_dir, exist_ok = True)
# ** save UMAP results
pd.DataFrame(ann.obsm['X_umap']).to_csv(
    os.path.join(l1_out_dir, "sa2_L1_umap.csv"),
    index = False, header = False
)
# ** save leiden results
ann.obs['leiden'].to_csv(
    os.path.join(l1_out_dir, "scRNA_sa2_annoy_k50_r0.8_leiden.csv"), index = True, header = True)



