import os
import sys

import numpy as np
import pandas as pd
import scipy

import snapatac2 as sa2
import scanpy as sc
import anndata as ad

from tmpPypkg.utils import reimport
from tmpPypkg import globalvar as gv

reimport("tmpPypkg.globalvar")

old_gmat_sa2ann = sa2.read(
  os.path.join("/projects/ps-renlab2/szu/projects",
               "CEMBA2",
               "19.snap2_integration",
               "src/main/resource",
               "sa2_gmat",
               "sa2_sa2default_ann_with_CPMlognorm.h5ad"),
  backed=None
)

ann_gmat = sc.AnnData(
  X=old_gmat_sa2ann.layers['sa2default_raw'],
  obs=old_gmat_sa2ann.obs
)
ann_gmat.var_names = old_gmat_sa2ann.var_names

ann_gmat.write_h5ad("tmp.snATAC.gmat.sa2ann.h5ad")

# * add snATAC meta
atacMeta = pd.read_csv(gv.ataccellmetafnm,
                       sep="\t",
                       header=0)
atacMeta.set_index("barcode2", drop=False, inplace=True)
old_obs = ann_gmat.obs
new_obs = atacMeta.loc[old_obs.index]
ann_gmat.obs = new_obs[['sample', 'barcode',
                   'MajorRegion', 'RegionName',
                   'subclass_id_label_v3',
                   'class_id_label_v3',
                   'ptgroup',
                   'NT_v3']].copy()
ann_gmat.write_h5ad(
  os.path.join(gv.pt_projd, "data/snATAC",
               "snATAC_gmat_obsv9.8_ann.h5ad")
)

# new_obs.loc[
#  new_obs.subclass_id_label_v3 == "326 OPC NN", ].gpL4.value_counts()


# check 4D and 4E
# OGC in 4E: 1368
# OGC in 4D: 1117
# OGC in 4H: 3162
sc_4H = new_obs.loc[
  new_obs.DissectionRegion == "4H", ].subclass_id_label_v3.value_counts()
sc_4D = new_obs.loc[
  new_obs.DissectionRegion == "4D", ].subclass_id_label_v3.value_counts()
sc_4E = new_obs.loc[
  new_obs.DissectionRegion == "4E", ].subclass_id_label_v3.value_counts()

# 70
sc_4H.index.intersection(sc_4D.index).shape
# 49
sc_4H.index.intersection(sc_4E.index).shape

sc_4D_4H = sc_4D[sc_4D.index.intersection(sc_4H.index)]
sc_4H_4D = sc_4H[sc_4D.index.intersection(sc_4H.index)]
# 0.68
scipy.stats.spearmanr(sc_4D_4H, sc_4H_4D)

sc_4E_4H = sc_4E[sc_4E.index.intersection(sc_4H.index)]
sc_4H_4E = sc_4H[sc_4E.index.intersection(sc_4H.index)]
# 0.29
scipy.stats.spearmanr(sc_4E_4H, sc_4H_4E)

# from bam file correlation
# 4D(CP-1) is closed to ACB
# 4E(ACB-2) is closed to CP

# So in snATAC: 4D should be 4E

sc_7H = new_obs.loc[
  new_obs.DissectionRegion == "7H", ].subclass_id_label_v3.value_counts()

sc_8H = new_obs.loc[
  new_obs.DissectionRegion == "8H", ].subclass_id_label_v3.value_counts()

sc_9G = new_obs.loc[
  new_obs.DissectionRegion == "9G", ].subclass_id_label_v3.value_counts()

# MEA related cells
sc_MEA = [i for i in sc_8H.index if "MEA" in i]
# same as above
# sc_MEA_7H = [i for i in sc_7H.index if "MEA" in i]
# 2977
n_MEA_7H = sc_7H[sc_MEA].sum()
# 4823
n_MEA_8H = sc_8H[sc_MEA].sum()
# 4035
n_MEA_9G = sc_9G[sc_MEA].sum()

# * limit snATAC to pt dissection
pt2disct = pd.read_csv(
  os.path.join(gv.pt_projd, "meta", "pairedtag_region2dissect.csv"),
  sep=",", header=None)

all_disct = ';'.join(
  pt2disct.iloc[:, 1].to_list()).split(";")

# temporally fit current dissect info in snATAC-seq
# so only 4D instead of 4E in the data
all_disct.remove("4E")
all_disct.append("4D")

all_disct = set(all_disct)

barcodes = new_obs.loc[
  new_obs.DissectionRegion.isin(all_disct), "barcode2"]
ann_gmat_pt = ann_gmat[barcodes, ].copy()

# * seperate NN and Neuron into two files
ann_gmat_pt_nn = ann_gmat_pt[ann_gmat_pt.obs.NT_v3 == "NN", ].copy()
ann_gmat_pt_neu = ann_gmat_pt[ann_gmat_pt.obs.NT_v3 != "NN", ].copy()

# save data
ann_gmat_pt.write_h5ad(
  os.path.join(gv.pt_projd, "data/snATAC",
               "snATAC_gmat_obsv9.8_ptregion_ann.h5ad")
)

ann_gmat_pt_nn.write_h5ad(
  os.path.join(gv.pt_projd, "data/snATAC",
               "snATAC_gmat_obsv9.8_ptregion_nn_ann.h5ad")
)

ann_gmat_pt_neu.write_h5ad(
  os.path.join(gv.pt_projd, "data/snATAC",
               "snATAC_gmat_obsv9.8_ptregion_neu_ann.h5ad")
)
