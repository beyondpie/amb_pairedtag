import os
from typing import List, Tuple
from multiprocessing import Pool
import numpy as np
import pandas as pd
import polars

import pyprojroot
import scanpy as sc
import anndata as ad
import snapatac2 as sa2

# * configs
proj_root = pyprojroot.here()
out_dir = os.path.join(proj_root, '01.clustering', "out")
scRNAseq_sa2_dir = os.path.join(out_dir, "scRNAseq_sa2")
# align genes: all the samples have the same genes
scRNAseq_sa2_aligene_dir = os.path.join(out_dir, "scRNAseq_sa2_aligene")
os.makedirs(scRNAseq_sa2_dir, exist_ok=True)
os.makedirs(scRNAseq_sa2_aligene_dir, exist_ok=True)
scRNAseq_sa2_anns = os.path.join(out_dir, "scRNAseq_sa2.anns.h5ad")
scRNAseq_sa2_all_ann = os.path.join(out_dir, "scRNAseq_sa2_all.ann.h5ad")

# * load cells' metadata
cell_meta = pd.read_csv(
    os.path.join(proj_root, "meta", "pairedtag.meta.v1.csv"),
    sep = ",",
    header = 0,
    index_col = None
)
cell_meta.set_index("barcode", inplace=True, drop=False)

# * save meta for one anndata
test_dir = os.path.join(proj_root, "01.clustering", "src", "test")
os.makedirs(test_dir, exist_ok = True)
i = 0
r = sa2.read(
    os.path.join(out_dir, "scRNAseq_sa2_aligene",
                 f"exp{i}_scRNA.aligene.sa2.h5ad"),
    backed = 'r')
rr = r.copy(filename = os.path.join(test_dir, "test1.h5ad"))
r.close()
rr.close()

r = sa2.read(
    os.path.join(test_dir, "test1.h5ad"),
    backed = 'r+'
)
# this will lead later loading r having AnnDataRead Error:
# AnnDataReadError: Above error raised while reading key '/' of type <class 'h5py._hl.files.File'> from /.

r_meta = polars.DataFrame(
    {"nCount_RNA": cell_meta.loc[r.obs_names, "nCount_RNA"].to_numpy(dtype = "u4")}
)
r.obsm["nCount_RNA"] = r_meta
r.close()

r = sa2.read(
    os.path.join(test_dir, "test1.h5ad"),
    backed = None
)


i = 0
r = sa2.read(
    os.path.join(out_dir, "scRNAseq_sa2_aligene",
                 f"exp{i}_scRNA.aligene.sa2.h5ad"),
    backed = 'r')
rr = r.copy(filename = os.path.join(test_dir, "test1.h5ad"))
r.close()
rr.close()

r = sa2.read(
    os.path.join(test_dir, "test1.h5ad"),
    backed = 'r+'
)
r_meta = polars.DataFrame(
    {"index": r.obs_names,
     "nCount_RNA": cell_meta.loc[r.obs_names, "nCount_RNA"].to_numpy(dtype = "u4")}
)
# unknown library error
r.obsm['index'] = polars.DataFrame({"index": r.obs_names})
