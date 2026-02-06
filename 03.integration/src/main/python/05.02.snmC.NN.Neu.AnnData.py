import os
import anndata as ad
from anndata.utils import make_index_unique
import numpy as np
import pandas as pd

import tmpPypkg.globalvar as gv
from tmpPypkg.singlecell import remove_black_list_region, remove_chromosomes
from tmpPypkg.singlecell import log_scale, neg_snmC

# from tmpPypkg.utils import show


# * meta
chrom_to_remove = ["chrM", "chrL"]
blv2fnm = os.path.join(gv.pt_projd, "meta", "mm10-blacklist.v2.bed")
snmCd = os.path.join(gv.pt_projd, "data/snmC_snm3C")

# * main
# load ann data
ann_mCG = ad.read_h5ad(
    filename=os.path.join(snmCd, "snmC.mCG.pt.Xonly.h5ad"), backed=None
)

ann_mCH = ad.read_h5ad(
    filename=os.path.join(snmCd, "snmC.mCH.pt.Xonly.h5ad"), backed=None
)

# load meta data
snmCmeta = pd.read_csv(
    os.path.join(snmCd, "raw", "CEMBA.mC.Metadata.csv"), sep=",", header=0
)
snmCmeta.set_index("cell", drop=False, inplace=True)

# transform the gene id to gene symbol
# geneNames = pd.read_csv(
#     os.path.join(gv.pt_projd, "meta", "gencode_GRCm38vM22_id2name.csv"),
#     sep=",",
#     header=0,
# )
# geneNames.set_index("gene_id", drop=False, inplace=True)
# ugnames = make_index_unique(geneNames.gene_name, join=".")
# geneNames.insert(2, "gene_name_uniq", ugnames)
# geneNames.to_csv(
#     os.path.join(gv.pt_projd, "meta",
#                  "gencode_GRCm38vM22_id2name2nameuniq.csv"),
#     sep=",",
#     index=False,
#     header=True,
# )

geneNames = pd.read_csv(
    os.path.join(gv.pt_projd, "meta", "gencode_GRCm38vM22_id2name2nameuniq.csv"),
    sep = ",",
    header = 0)
geneNames.set_index("gene_id", drop=False, inplace=True)


old_varnames = ann_mCG.var_names
varnames = geneNames.loc[old_varnames].gene_name_uniq
ann_mCG.var_names = varnames
ann_mCH.var_names = varnames

# add meta data to ann
ann_mCG.obs = snmCmeta.loc[ann_mCG.obs_names]
ann_mCH.obs = snmCmeta.loc[ann_mCH.obs_names]

ann_mCG.write_h5ad(os.path.join(snmCd, "snmC.mCG.pt.h5ad"))
ann_mCH.write_h5ad(os.path.join(snmCd, "snmC.mCH.pt.h5ad"))

# mix mCG and mCH to provide the approximate gene mat
classes = snmCmeta.SubClass.apply(lambda x: x.split(" ")[-1])
# 23966 cells
cellNN = snmCmeta.loc[classes.isin(["NN"])].cell
# 1747 cells
cellIMN = snmCmeta.loc[classes.isin(["IMN"])].cell
# 7669 cells
cellDGlut = snmCmeta[snmCmeta.SubClass == "DG Glut"].cell
cell_NN_IMN_DGlut = pd.concat([cellNN, cellIMN, cellDGlut])

cell_usemCG = ann_mCG.obs_names[ann_mCG.obs_names.isin(cell_NN_IMN_DGlut)]
cell_usemCH = ann_mCH.obs_names[~ann_mCG.obs_names.isin(cell_usemCG)]

ann_mCG_NNIDG = ann_mCG[cell_usemCG, :]
ann_mCH_Neu = ann_mCH[cell_usemCH, :]

ann_snmCgmat = ad.concat([ann_mCG_NNIDG, ann_mCH_Neu])
ann_snmCgmat.var = ann_mCG.var

ann_snmCgmat.obs.insert(22, "isNN", False)
# 5024 NN
ann_snmCgmat.obs.isNN[ann_snmCgmat.obs_names.isin(cellNN)] = True
ann_snmCgmat.write_h5ad(os.path.join(snmCd, "snmC.pt.mCH_mCG.for.gmat.h5ad"))

# downsample the classes based on the clustering results
ann_snmCgmat = ad.read_h5ad(
    os.path.join(snmCd, "snmC.pt.mCH_mCG.for.gmat.h5ad"))
neu_ann = ann_snmCgmat[~ann_snmCgmat.obs.isNN,]

nds = 1500
neu_dscell = (
    neu_ann.obs.groupby("SubClass", observed=True)
    .apply(lambda x: x.sample(min(nds, x.shape[0])))
    .cell
)
nds = 500
nn_ann = ann_snmCgmat[ann_snmCgmat.obs.isNN,]
nn_dscell = nn_ann.obs.groupby("SubClass", observed=True).apply(
    lambda x: x.sample(min(nds, x.shape[0])).cell
)

# and limit to k8
k8g = pd.read_csv(
    os.path.join(gv.pt_projd, "meta", "AIT21_k8_markers.txt"), header=None
).iloc[:, 0]
ann_k8g = neu_ann.var_names.isin(k8g)
neu_dsann = neu_ann[neu_dscell, ann_k8g].copy()
nn_dsann = nn_ann[nn_dscell, ann_k8g].copy()

# check mC or neg mC
# g = ["Snap25", "Gad1", "Gad2", "Atp1a2", "Mbp"]
# scs = ['Oligo NN', 'Astro-TE NN', 'CA3 Glut', 'DG Glut',
#      'DG-PIR Ex IMN', 'Lamp5 Gaba']
# x = ann_snmCgmat[:, g]
# a = x.obs.groupby("SubClass", observed=True).apply(
#     lambda i: np.mean(x[i.cell, :].X, axis=0))
# g_snap25 = pd.Series([i[0] for i in a[scs]], index=scs)
# g_Gad1 = pd.Series([i[1] for i in a[scs]], index=scs)
# g_Gad2 = pd.Series([i[2] for i in a[scs]], index=scs)
# g_Atp1a2 = pd.Series([i[3] for i in a[scs]], index=scs)
# g_Mbp = pd.Series([i[4] for i in a[scs]], index=scs)

# filter chromesome and blacklist
remove_chromosomes(nn_dsann, exclude_chromosomes=chrom_to_remove)
remove_chromosomes(neu_dsann, exclude_chromosomes=chrom_to_remove)
remove_black_list_region(nn_dsann, black_list_path=blv2fnm)
remove_black_list_region(neu_dsann, black_list_path=blv2fnm)

# log_scale and take negative value then
log_scale(nn_dsann, with_mean=True)
log_scale(neu_dsann, with_mean=True)
neg_snmC(nn_dsann)
neg_snmC(neu_dsann)

# save anndata for integration
nn_dsann.write_h5ad(os.path.join(snmCd, "snmC.scaled.gmat.nn.ds.h5ad"))
neu_dsann.write_h5ad(os.path.join(snmCd, "snmC.scaled.gmat.neu.ds.h5ad"))


# * 2024-11-13
# use all the data without downsampling
k8neuann = neu_ann[ :, ann_k8g].copy()
k8nnann = nn_ann[ :, ann_k8g].copy()

# filter chromesome and blacklist
remove_chromosomes(k8nnann, exclude_chromosomes=chrom_to_remove)
remove_chromosomes(k8neuann, exclude_chromosomes=chrom_to_remove)
remove_black_list_region(k8nnann, black_list_path=blv2fnm)
remove_black_list_region(k8neuann, black_list_path=blv2fnm)

# log_scale and take negative value then
log_scale(k8nnann, with_mean=True)
log_scale(k8neuann, with_mean=True)
neg_snmC(k8nnann)
neg_snmC(k8neuann)

# save anndata for integration
k8nnann.write_h5ad(os.path.join(snmCd, "snmC.scaled.gmat.k8.nn.h5ad"))
k8neuann.write_h5ad(os.path.join(snmCd, "snmC.scaled.gmat.k8.neu.h5ad"))

# * 2024-11-18
# merge nn and neu
ann_nnds = ad.read_h5ad(os.path.join(snmCd, "snmC.scaled.gmat.nn.ds.h5ad"),
                        backed=None)
ann_neuds = ad.read_h5ad(os.path.join(snmCd, "snmC.scaled.gmat.neu.ds.h5ad"),
                         backed=None)
ann_ds = ad.concat([ann_nnds, ann_neuds], merge="same")
ann_ds.write_h5ad(os.path.join(snmCd, "snmC.scaled.gmat.k8.ds.h5ad"))

