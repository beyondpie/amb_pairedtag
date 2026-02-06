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

# * functions
def sa2_read10X(mtx_dir:str, outfnm:str) -> sa2.AnnData:
    """
    mtx_dir: mtx_dir
    outfnm: output file name
    """
    import scanpy as sc
    import snapatac2 as sa2
    if not os.path.exists(mtx_dir):
        raise ValueError(f"mtx_dir does not exist: {mtx_dir}")
    if os.path.exists(outfnm):
        raise ValueError(f"outfnm already exists: {outfnm}")
    t = sc.read_10x_mtx(path = mtx_dir)
    r = sa2.AnnData(filename = outfnm)
    r.X = t.X
    r.obs_names = t.obs_names
    r.var_names = t.var_names
    return r

# * read csv file
meta = pd.read_csv(
    '../resource/scRNAseq_samples.passQC.231120.csv',
    sep = ",",
    header = 0,
    index_col=None)

# * mtx to sa2 AnnData
for i in range(meta.shape[0]):
    r = sa2_read10X(
        mtx_dir = meta.loc[i, "mtx_dir"],
        outfnm = os.path.join(scRNAseq_sa2_dir, f"exp{i}_scRNA.sa2.h5ad"))
    r.close()
    print(f"exp{i} done")

# * get genes across all samples
r = sa2.read(
    filename = os.path.join(scRNAseq_sa2_dir, "exp0_scRNA.sa2.h5ad"),
    backed = "r")
genes: List[str] = r.var_names
r.close()
for i in range(1, meta.shape[0]):
    print(f"loading exp{i}.")
    r = sa2.read(
        filename = os.path.join(scRNAseq_sa2_dir, f"exp{i}_scRNA.sa2.h5ad"),
        backed = "r")
    genes = list(set(genes) & set(r.var_names))
    r.close()

# * load cells' metadata
cell_meta = pd.read_csv(
    os.path.join(proj_root, "meta", "pairedtag.meta.v1.csv"),
    sep = ",",
    header = 0,
    index_col = None
)
cell_meta.set_index("barcode", inplace=True, drop=False)
barcodes_set = set(cell_meta.barcode)

# * subset AnnData based on genes and cells
def reduce_expi(i: int) -> None:
    print(f"load exp{i}.")
    r1 = sa2.read(
        filename = os.path.join(scRNAseq_sa2_dir, f"exp{i}_scRNA.sa2.h5ad"),
        backed = "r")
    barcodes_r2 = [i for i in r1.obs_names if i in barcodes_set]
    r2 = r1.subset(
        obs_indices = barcodes_r2,
        var_indices = genes,
        out = os.path.join(scRNAseq_sa2_aligene_dir,
                           f"exp{i}_scRNA.aligene.sa2.h5ad"))
    r1.close()
    r2.close()
    print(f"reduce exp{i} done.")

with Pool(3) as pool:
    pool.map(reduce_expi, range(meta.shape[0]))


# * merge all AnnData to AnnDataSet
files:List[Tuple[str, str]] = [
    (f"exp{i}", os.path.join(scRNAseq_sa2_aligene_dir,
                             f"exp{i}_scRNA.aligene.sa2.h5ad"))
         for i in range(meta.shape[0])]
anns = sa2.AnnDataSet(adatas = files, filename = scRNAseq_sa2_anns)
anns.close()

# * add meta data
# FIXME: we can add meta in the following way,
# but we then cannot load it to memory (load it with the backed 'r' is fine).
# So we now just put meta data in a seperated file.

# anns = sa2.read_dataset(filename = scRNAseq_sa2_anns,
#                         update_data_locations = None,
#                         mode = 'r+')

# anns_meta = polars.DataFrame(
#     {
#         "nCount_RNA": cell_meta.loc[
#             anns.obs_names, "nCount_RNA"].to_numpy(dtype = "u4"),
#         "nFeature_RNA": cell_meta.loc[
#             anns.obs_names, "nFeature_RNA"].to_numpy(dtype = "u4"),
#         "percent.mt": cell_meta.loc[
#             anns.obs_names, "percent.mt"].to_numpy(dtype = "f4"),
#         "brainregion": cell_meta.loc[
#             anns.obs_names, "brainregion"].to_numpy(dtype = "U32"),
#         "modality": cell_meta.loc[
#             anns.obs_names, "modality"].to_numpy(dtype = "U32"),
#         "oriBarcode": cell_meta.loc[
#             anns.obs_names, "oriBarcode"].to_numpy(dtype = "u4"),
#         "sublib": cell_meta.loc[
#             anns.obs_names, "sublib"].to_numpy(dtype = "U32"),
#         "sex": cell_meta.loc[anns.obs_names, "sex"].to_numpy(dtype = "U8"),
#         "rep": cell_meta.loc[anns.obs_names, "rep"].to_numpy(dtype = "U8"),
#         "exp": cell_meta.loc[anns.obs_names, "exp"].to_numpy(dtype = "U32")
#     }
# )
# for i in anns_meta.columns:
#     anns.obsm[i] = polars.DataFrame({i : anns_meta[i]})
# anns.close()

# * AnnDataSet to one AnnData
anns = sa2.read_dataset(filename = scRNAseq_sa2_anns,
                        mode = 'r')
r = anns.to_adata(copy_x = True, file = scRNAseq_sa2_all_ann)
anns.close()
r.close()

# # * check meta data
# orig_meta = cell_meta
# ann_all = sa2.read(filename = scRNAseq_sa2_all_ann, backed = 'r')
# meta_ann_all = ann_all.obsm[
#     ['rep', 'sex', 'sublib', 'exp', 'brainregion',
#      'modality', 'percent.mt', 'oriBarcode',
#      'nFeature_RNA', 'nCount_RNA']].to_pandas()


# acol = 'nFeature_RNA'

# a = meta_ann_all[acol].to_pandas().squeeze()
# a.index = ann_all.obs_names
# b = cell_meta.loc[a.index, acol]

# a = a.astype('int64')
# a.equals(cell_meta.loc[a.index, acol])

# * prepare test data for clustering
import random
test_ann_outfnm = os.path.join(proj_root, "01.clustering",
                      "src/test/resource",
                      "scRNAseq_sa2_test.ann.h5ad")
os.makedirs(os.path.dirname(test_ann_outfnm), exist_ok = True)
ann = sa2.read(filename = scRNAseq_sa2_all_ann, backed = 'r')
rdm_index = random.sample(population = range(ann.shape[0]), k = 20000)
sub_ann = ann.subset(
    obs_indices = rdm_index, out = test_ann_outfnm)
barcodes:List[str] = sub_ann.obs_names
sub_ann.close()

L2: List[int] = [0] * 10000 + [1] * 10000
b2g: pd.DataFrame = pd.DataFrame(
    data = {'barcode': barcodes,
            'L1': [0]*len(barcodes),
            'L2': L2})
b2g.to_csv(
    os.path.join(os.path.dirname(test_ann_outfnm), "b2g_test.csv"),
    sep = ",",
    header = True,
    index = False
)

# * prepare L1 clustering meta
ann = sa2.read(filename = scRNAseq_sa2_all_ann, backed = 'r')
barcodes: pd.DataFrame = pd.DataFrame(data = ann.obs_names,
                        index = ann.obs_names,
                        columns = ['barcode'])
barcodes['L0'] = '0'
barcodes.to_csv(
    os.path.join(proj_root,
                 "01.clustering",
                 "src/main/resource",
                 "barcode2group_L0.csv"),
    sep = ",",
    header = True, index = False
)

# * re-prepare the h5ad and meta file for L1 clustering
cell_meta = pd.read_csv(
    os.path.join(proj_root, "meta", "pt.barcode.meta.with.qc.v2.csv"),
    sep = ",",
    header = 0
)

barcode: pd.DataFrame = cell_meta[
    ~cell_meta.isLQ].barcode.to_frame(name = "barcode")
barcode.insert(1, "L0", 0)

barcode.to_csv(
    os.path.join(proj_root, "01.clustering", "src/main/resource",
                 "barcode2group_L0_afqc.csv"),
    sep = ",",
    header = True, index = False
)

# prepare anndata in the meanwhile
ann: sa2.AnnData = sa2.read(
    filename = os.path.join(proj_root, "01.clustering",
                            "out", "scRNAseq_sa2_all.ann.h5ad"),
    backed = 'r')

ann_qc = ann.subset(obs_indices=barcode.barcode,
                    out = os.path.join(proj_root, "01.clustering",
                                       "out",
                                       "scRNAseq_sa2_afqc.ann.h5ad"))

ann_qc = sa2.read(
    filename = os.path.join(proj_root, "01.clustering",
                            "out",
                            "scRNAseq_sa2_afqc.ann.h5ad"),
    backed = 'r'
)
ann_qc.close()
