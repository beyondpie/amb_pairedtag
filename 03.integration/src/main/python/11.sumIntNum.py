"""
Summarize number of cells used for integration.
"""

import sys
import os
import re
import snapatac2 as sa2
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import tmpPypkg.globalvar as gv

# * configs
projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"

# * main

# get number of cells per region for different modalities


# 1. ATAC-seq
# See 06.snATAC.NN.Neu.AnnData.py for details
snATACd = os.path.join(gv.pt_projd, "data/snATAC")
ann_nn = ad.read_h5ad(
    os.path.join(snATACd, "snATAC_gmat_obsv9.8_pteregion_nn_ann.h5ad"),
    backed='r')
ann_neu = ad.read_h5ad(
    os.path.join(snATACd, "snATAC_gmat_obsv9.8_ptregion_neu_ann.h5ad"),
    backed='r'
)

# 472528 cells
nATAC = ann_nn.shape[0]  + ann_neu.shape[0]

# 2. AllenRNA-seq
with open(
    os.path.join(projd, "data", "allen",
                 "allen_10xv3_region2count.csv"), 'r') as f:
    lines = [l.strip() for l in f.readlines()]
    allen_10xv3_r2cnt = {
       l.split(",")[0]: int(l.split(",")[1]) for l in lines
    }
# based on resources/Allen_reference_data_locator_ZW_20240131.csv
pt2allen_region = {
    'CPU': ['STR - STRd'],
    'HYP': ['HY -  MEZ-PVZ-PVR', 'HY LZ', 'HY - HYmm', 'HY - HYml',
            'HY - HYm2', 'HY - HYa1', 'HY - HYa2', 'HY - HYpm', 'HY - HYpl'],
    'HC': ['HIP - CA', 'HIP'],
    'ENT': ['ENT'],
    'AMY': ['STR - sAMY'],
    'NAC': ['STR - STRv'],
    'VTA_SnR': ['MB - VTA-SN'],
    'PFC': ['ACA', 'PL-ILA-ORB']
}

nallen = 0
for k, v in pt2allen_region.items():
    for r in v:
        nallen += allen_10xv3_r2cnt[r]
# nallen: 656346
    
# 3. snmC
# see 05.02.snmC.NN.Neu.AnnData.py for details
snmCd = os.path.join(gv.pt_projd, "data/snmC_snm3C")
ann_snmCgmat = ad.read_h5ad(
    os.path.join(snmCd, "snmC.pt.mCH_mCG.for.gmat.h5ad"), backed = 'r' )
# 63792 cells

# * downsampel for co-emebedding
# check 01.split.Allen.py last lines for details.
