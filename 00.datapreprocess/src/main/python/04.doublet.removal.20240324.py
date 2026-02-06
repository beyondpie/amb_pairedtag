import sys
import os
from multiprocessing import Pool
from typing import List, Tuple, Dict
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import snapatac2 as sa2

import pyprojroot

proj_root = pyprojroot.here()

# * load estimation of doublet rates.
dltrate_1to13: List[float] = [
    0.377,
    0.232,
    0.485,
    0.133,
    0.898,
    5.234,
    0.551,
    2.323,
    1.330,
    1.314,
    0.762,
    3.247,
    0.933,
]
dltrate_14: float = round(np.mean(dltrate_1to13), 3)
dltrate_all = pd.DataFrame(
    np.around(np.concatenate([dltrate_1to13, [dltrate_14]]) * 0.01, 5),
    columns=["dltr"]
)
dltrate_all.insert(loc=0,
                   column="exp",
                   value=[f"exp{i}" for i in range(1, 15)])
dltrate_all.to_csv(
    os.path.join(proj_root, "meta", "est.dlt.rate.20240324.csv"),
    sep=",",
    header=True,
    index=False
)

# * set up AnnData
mtx_dir = os.path.join("/projects/ps-renlab2/zhw063",
                       "99.MouseBrainPairedTag",
                       "2024Mar_new_mtxs")


def load_RNA_mtx(i: int) -> ad.AnnData:
    path = os.path.join(mtx_dir, f"AMB_Exp{i}_RNA")
    print(f"read RNA mtx from {path}.")
    return sc.read_10x_mtx(path=path)


anns: List[ad.AnnData] = [load_RNA_mtx(i) for i in range(1, 15)]
adatas: Dict[str, ad.AnnData] = {
    f"Exp{i+1}": a for i, a in enumerate(anns)
}
ann_all = ad.concat(adatas, label="Exp")
ann_all.write_h5ad(
    filename=os.path.join(proj_root, "data",
                          "pairedtag_ann",
                          "pt.rawRNAseq.20240324.h5ad"))

# * doublet score/prob from SnapATAC2
# load the raw snRNA-seq part
ann_all: ad.AnnData = ad.read_h5ad(
    filename=os.path.join(proj_root, "data",
                          "pairedtag_ann",
                          "pt.rawRNAseq.20240324.h5ad"),
    backed=None)
# load cell meta (cells passing mito filtering)
cellmeta: pd.DataFrame = pd.read_csv(
    os.path.join(proj_root, "02.track",
                 "src/main/resource",
                 "cellMeta.20240321.csv"),
    sep=",",
    header=0
)
# filter ann_all by cellmeta
ann: ad.AnnData = ann_all[cellmeta.barcode, :]
ann.write_h5ad(
    filename=os.path.join(proj_root, "data",
                          "pairedtag_ann",
                          "pt.RNAseq.filterbymito.20240327.h5ad")
)
del ann_all
del ann

# load doublet rate per exp
dlts: pd.DataFrame = pd.read_csv(
    os.path.join(proj_root, "meta", "est.dlt.rate.20240324.csv"),
    sep=",", header=0
)
dlts.set_index("exp", drop=False, inplace=True)

# get snapatac2 anns
mtx_dir = os.path.join("/projects/ps-renlab2/zhw063",
                       "99.MouseBrainPairedTag",
                       "2024Mar_new_mtxs")
sa2ann_dir = os.path.join(proj_root, "data",
                          "pairedtag_ann",
                          "snapatac2")


def load_RNAmtx_sa2ann(i: int) -> None:
    path = os.path.join(mtx_dir, f"AMB_Exp{i}_RNA")
    topath = os.path.join(sa2ann_dir, f"sa2.AMB_Exp{i}_RNA.h5ad")
    print(f"read RNA mtx from {path} to {topath}.")
    if os.path.exists(topath):
        print(f"{topath} exists, and remove it.")
        os.remove(path=topath)
    r = sa2.read_10x_mtx(path=path,
                         file=os.path.join(sa2ann_dir,
                                           f"sa2.AMB_Exp{i}_RNA.h5ad"))
    # genes: pd.DataFrame = pd.read_csv(
    #     os.path.join(path, "genes.tsv"),
    #     header = None, sep = "\t"
    # )
    # var_names = ad.utils.make_index_unique(
    #     index = pd.Index(genes[1].values)).to_list()
    # r.var_names = var_names
    r.close()
    return None


for i in range(1, 15):
    load_RNAmtx_sa2ann(i)

anns: List[Tuple[str, sa2.AnnData]] = [
    (f"exp{i}",
     sa2.read(os.path.join(sa2ann_dir, f"sa2.AMB_Exp{i}_RNA.h5ad"),
              backed='r')
     )
    for i in range(1, 15)
]

# get unified genes
list_of_varnames = [r.var_names for _, r in anns]
for _, a in anns:
    a.close()

varnames = list_of_varnames[0]
for i, v in enumerate(list_of_varnames):
    print(f"{i}")
    varnames = [g for g in v if g in varnames]

# get anns with unified var_names


def unify_ann(i: int) -> None:
    f1 = os.path.join(sa2ann_dir, f"sa2.AMB_Exp{i}_RNA.h5ad")
    print(f"Unifying ann {f1}.")
    f2 = os.path.join(sa2ann_dir,
                      f"sa2.AMB_Exp{i}_RNA.unified.genesymol.h5ad")
    if os.path.exists(f2):
        os.remove(f2)
    r1: sa2.AnnData = sa2.read(f1, backed='r')
    r2 = r1.subset(var_indices=varnames, out=f2)
    r2.close()
    return None


for i in range(1, 15):
    unify_ann(i)

# still ENSMUG, but they are aligned to the same dimension
annfnms: List[Tuple[str, sa2.AnnData]] = [
    (f"exp{i}",
     sa2.read(os.path.join(proj_root, "data/pairedtag_ann",
                           "snapatac2",
                           f"sa2.AMB_Exp{i}_RNA.unified.genesymol.h5ad"),
              backed='r'))
    for i in range(1, 15)
]

annset = sa2.AnnDataSet(
    adatas=annfnms, add_key='exp',
    filename=os.path.join(proj_root, "data/pairedtag_ann",
                          "pt.rawRNAseq.sa2annset.h5ad")
)

# AMB_Exp1_RNA genes.tsv has some bug:
# some ENSUMG are not unique.
genes: pd.DataFrame = pd.read_csv(
    os.path.join("/projects/ps-renlab2/zhw063",
                 "99.MouseBrainPairedTag",
                 "2024Mar_new_mtxs", "AMB_Exp2_RNA", "genes.tsv"),
    sep="\t", header=None)
genes.columns = ["en", "sym"]
genes.drop_duplicates(keep="first", inplace=True)
genes.set_index("en", drop=False, inplace=True)
genes = genes.loc[varnames]
gene_symbols = ad.utils.make_index_unique(
    pd.Index(genes.sym.values)).to_list()
genes['gene_symbol'] = gene_symbols

genes.to_csv(os.path.join(proj_root, "meta", "Mouse.ENSMUG2genesymbol.tsv"),
             sep="\t",
             header=True, index=False)
g = genes.loc[annset.var_names].gene_symbol.to_list()
annset.var_names = g

ann = annset.to_adata(
    file=os.path.join(
        proj_root, "data/pairedtag_ann",
        "pt.rawRNAseq.sa2ann.genesymbol.h5ad"))

# load cell meta (cells passing mito filtering)
cellmeta: pd.DataFrame = pd.read_csv(
    os.path.join(proj_root, "02.track",
                 "src/main/resource",
                 "cellMeta.20240321.csv"),
    sep=",",
    header=0
)
# filtered by mito
ann = sa2.read(
    os.path.join(proj_root, "data/pairedtag_ann",
                 "pt.rawRNAseq.sa2ann.genesymbol.h5ad"),
    backed='r'
)
ann_fm = ann.subset(obs_indices=cellmeta.barcode.to_list(),
                    out=os.path.join(proj_root, "data/pairedtag_ann",
                                     "pt.RNAseq.sa2ann.filterbymito.h5ad"))
ann.close()
ann_fm.close()

# * get individiual ann after filtering cells.
sa2ann_dir = os.path.join(proj_root, "data/pairedtag_ann",
                          "snapatac2")
def filter_cell(i: int) -> None:
    infnm = os.path.join(sa2ann_dir,
                         f"sa2.AMB_Exp{i}_RNA.unified.genesymol.h5ad")
    # filter cells based on cell meta
    outfnm = os.path.join(sa2ann_dir,
                          f"sa2.AMB_Exp{i}_RNA.bfdlt.h5ad")
    print(f"from {infnm} to {outfnm}.")
    if os.path.exists(outfnm):
        os.remove(outfnm)
    r1: sa2.AnnData = sa2.read(infnm, backed='r')
    r2 = r1.subset(
        obs_indices=pd.Series(r1.obs_names).isin(cellmeta.barcode).to_list(),
        out=outfnm
    )
    r2.close()
    r1.close()
    return None


for i in range(1, 15):
    filter_cell(i)


