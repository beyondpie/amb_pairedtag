import sys
import os
from typing import List, Tuple
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import snapatac2 as sa2

import pyprojroot

proj_root = pyprojroot.here()


# * load doublet rate for different experiments.
exp2dltr: pd.DataFrame = pd.read_csv(
    os.path.join(proj_root, "meta", "doublet.rate.est.csv"),
    sep=",",
    header=0,
)
exp2dltr.set_index("exp", drop=False, inplace=True)

# * load snapatac2 object
# exp0 - exp13 in obs['sample']
sa2ann_rna = sa2.read(
    os.path.join(
        proj_root, "01.clustering",
        "out", "nodlt", "L1", "RNA_k8_npc50_k50_L0_0.h5ad"
    ),
    backed="r")

# * performa doublet rate estimation per sample
ann_inmem: ad.AnnData = sa2ann_rna.to_memory()


def get_dlt(annexp_i: int) -> pd.DataFrame:
    annexp = f"exp{annexp_i}"
    ptexp = f"Exp{annexp_i + 1}"
    dltr = exp2dltr.loc[ptexp, "value"]
    print(f"run scrublet on ann exp: exp{annexp_i}.")
    print(f"init doublet rate is: {dltr}.")
    sub_ann = ann_inmem[ann_inmem.obs["sample"].isin([annexp])]
    dlt_result: Tuple[np.ndarray, np.ndarray] = sa2.preprocessing.scrublet(
        sub_ann, features=None, expected_doublet_rate=dltr, inplace=False
    )
    r = pd.DataFrame(
        data=np.transpose(np.vstack(dlt_result)),
        index=sub_ann.obs_names,
        columns=["dlt_prob", "dlt_score"],
    )
    return r


# validate the order of dlt_result
# b = ann_inmem[ann_inmem.obs['sample'].isin(['exp0'])].copy()
# sa2.preprocessing.scrublet(
#     b, features = None,
#     expected_doublet_rate=exp2dltr.loc['Exp1', 'value'],
#     inplace = True
# )
dlt_results: List[pd.DataFrame] = list(map(get_dlt, range(14)))
dlt_all: pd.DataFrame = pd.concat(dlt_results)
dlt_all.insert(0, "barcode", dlt_all.index, allow_duplicates=False)
dlt_all.to_csv(
    os.path.join(proj_root, "meta", "sa2_scrublet_8k_all.csv"),
    sep=",",
    header=True,
    index=False,
)

# * filter doublets
cellmeta: pd.DataFrame = pd.read_csv(
    os.path.join(proj_root, "meta", "pairedtag.meta.v1.csv"),
    sep = ",",
    header = 0,
    index_col = None
)
cellmeta.set_index("barcode", drop=False, inplace=True)
dlt_all: pd.DataFrame = pd.read_csv(
    os.path.join(proj_root, "meta", "sa2_scrublet_8k_all.csv"),
    sep=",",
    header=0,
    index_col=None)
dlt_all.set_index("barcode", drop=False, inplace=True)
dlt_all.insert(loc = 3, column = "exp",
               value = cellmeta.loc[dlt_all.index, "exp"])

# ** get exp-specific doublet prob cutoff
def get_dltprob(rows: pd.DataFrame) -> float:
    exp: str = rows.iloc[0]['exp']
    threshold: float = min(0.5,
                           rows.dlt_prob.quantile(
                               1 - exp2dltr.loc[exp].value))
    return threshold

def get_dlt_thres(
        rows: pd.DataFrame, oncol: str = 'dlt_prob') -> float:
    exp: str = rows.iloc[0]['exp']
    threshold: float = rows[oncol].quantile(
        1 - exp2dltr.loc[exp].value)
    return threshold

# BUG: problem at Exp8, threshold is 0.0
# Only about 10 cells has prob > 0.9, all the others < 0..002
# expecially lots of them < 10^(-17)
exp_dltprob_thres = dlt_all.groupby(by = 'exp').apply(
    func = get_dlt_thres, oncol = 'dlt_prob')
exp_dltscore_thres = dlt_all.groupby(by = 'exp').apply(
    func = get_dlt_thres, oncol = 'dlt_score')

exp_dlt_thres = pd.concat([exp_dltprob_thres, exp_dltscore_thres],
                          axis = 1)
exp_dlt_thres.columns = ["dltprob_thres", "dltscore_thres"]
exp_dlt_thres.insert(0, 'exp', exp_dlt_thres.index)
exp_dlt_thres.to_csv(
    os.path.join(proj_root, "meta",
                 "exp2dlthres.L1.k8.csv"),
    sep = ",",
    header = True,
    index = False
)
# now filter by dltscore
dlt_filtered = dlt_all.groupby(by = 'exp').apply(
    lambda x: x[x.dlt_score <= exp_dltscore_thres.loc[x.iloc[0].exp]]
)
dlt_filtered.to_csv(
    os.path.join(proj_root, "meta",
                 "barcode2dlt.filtered.L1.k8.dltscore.csv"),
    sep = ",",
    header = True,
    index = False
)
