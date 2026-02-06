import sys
import os
from typing import List
import math
import pandas as pd
import scanpy as sc
import anndata as ad
import snapatac2 as sa2

import pyprojroot

proj_root = pyprojroot.here()
exp2dltr: pd.DataFrame = pd.read_csv(
    os.path.join(proj_root, "meta", "est.dlt.rate.20240324.csv"),
    sep=",", header=0
)
exp2dltr.set_index("exp", drop=False, inplace=True)

dlt_dir = os.path.join(proj_root, "00.datapreprocess",
                       "out/sa2dlt")
# load doublets


def add_exp(i) -> pd.DataFrame:
    r: pd.DataFrame = pd.read_csv(os.path.join(dlt_dir,
                                               f"sa2.AMB_Exp{i}_RNA.dlt.csv"),
                                  sep=",",
                                  header=0)
    r.insert(0, "exp", value=f"exp{i}")
    return (r)


dlts: List[pd.DataFrame] = [add_exp(i) for i in range(1, 15)]

dlt: pd.DataFrame = pd.concat(dlts, axis=0)

# load cells used to define doublets
mtx_dir = os.path.join("/tscc/projects/ps-renlab2/zhw063",
                       "99.MouseBrainPairedTag",
                       "2024Mar_new_mtxs")
dlt_zhw: pd.DataFrame = pd.read_csv(
    os.path.join(mtx_dir, "20240325_brain.rna.metadata.txt"),
    header=0,
    sep=" "
)
dlt_zhw.insert(0, "barcode", dlt_zhw.index)

barcode_hm: List[str] = dlt_zhw[dlt_zhw.anno == "human"].index.to_list()
barcode_mix: List[str] = dlt_zhw[dlt_zhw.anno == "mix"].index.to_list()

# set dlt cutoff
outdir = os.path.join(proj_root, "meta")


def get_dlt_thres(
        rows: pd.DataFrame, oncol: str = 'dlt_prob') -> float:
    exp: str = rows.iloc[0]['exp']
    threshold: float = rows[oncol].quantile(
        1 - exp2dltr.loc[exp].dltr)
    return threshold


ub_probs = dlt.groupby(by='exp').apply(
    func=get_dlt_thres, oncol="doublet_probability"
)
ub_score: pd.Series = dlt.groupby(by='exp').apply(
    func=get_dlt_thres, oncol="doublet_score"
)


cell_afdlt: pd.DataFrame = dlt.groupby(by='exp').apply(
    lambda x: x[x.doublet_score <= ub_score.loc[x.iloc[0].exp]]
)
# save threshold
a1 = pd.DataFrame(ub_probs, columns=["sa2_top3k_scrublet_prob"])
a1.insert(0, "exp", ub_probs.index)
a1.insert(2, "sa2_top3k_scrublet_score", ub_score.loc[ub_probs.index])
a1.insert(3, "est_dlt", exp2dltr.loc[ub_probs.index].dltr)
a1.to_csv(os.path.join(outdir, "doublet.upbound.sa2.top3k.scrubelt.csv"),
          sep=",", header=True, index=False)

# generate cell meta with dlt info
# - 2 meta: one is raw, one is after filtering
# 1. raw cell meta
raw_cellmeta = pd.merge(dlt, dlt_zhw, on="barcode")
raw_cellmeta.insert(raw_cellmeta.shape[1], "isdlt",
                    ~ raw_cellmeta.barcode.isin(cell_afdlt.barcode))
raw_cellmeta.to_csv(
    os.path.join(outdir, "pt.barcode.meta.with.dlt.20240330.csv"),
    sep=",",
    index=False,
    header=True
)
# 2. get cell meta after filtering
filtered_cellmeta = raw_cellmeta[
    (raw_cellmeta.anno == 'mouse') &
    (~raw_cellmeta.isdlt)
]
filtered_cellmeta.to_csv(
    os.path.join(outdir, "pt.barcode.meta.after.dlt.20240330.csv"),
    sep=",",
    index=False,
    header=True
)

# generate sa2ann without dlt for L1 clustering
raw_ann: sa2.AnnData = sa2.read(
    os.path.join(proj_root, "data/pairedtag_ann",
                 "pt.RNAseq.sa2ann.filterbymito.h5ad"),
    backed = 'r'
)

obs_index: List[bool] = pd.Series(raw_ann.obs_names).isin(
    filtered_cellmeta.barcode).to_list()
ann4L1 = raw_ann.subset(obs_indices = obs_index,
                        out = os.path.join(proj_root, "data/pairedtag_ann",
                                           "pt.RNAseq.sa2ann4L1.h5ad"))
raw_ann.close()
ann4L1.close()

# prepare L1 meta
a = pd.DataFrame(filtered_cellmeta.barcode, columns = ["barcode"])
a.insert(1, "L0", 0)
a.to_csv(
    os.path.join(proj_root, "01.clustering",
                 "src/main/resource",
                 "pt.barcode2group.L0.20240330.csv"),
    header = True,
    index = False
)
