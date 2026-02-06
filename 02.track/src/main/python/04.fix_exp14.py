"""
For Issue "Correct cell meta bug on Experiment 14 #5"

For Exp14, barcodes #23 and #24 are mistakenly assigned to male PFC H3K9me3.
- #23 should be ERC-Male-RepA-H3K9me3
- #24 should be AMY-Male-RepA-H3K9me3.
"""

import os
import sys
import re

import pandas as pd

import pyprojroot

proj_dir = pyprojroot.here()
cell_meta_file = os.path.join(proj_dir, "meta",
                              "pairedtag.meta.v1.csv")


# * fix pairedtag.meta.v1.csv
cell_meta = pd.read_csv(cell_meta_file, sep = ",", header = 0,
                        index_col = False)
with open(cell_meta_file, 'r') as f:
    metalines = f.readlines()
    
pt_e14_b23 = "Exp14,\w\w\w:\w\w:\w\w:23"
row_index = [i-1
    for i,l in enumerate(metalines)
    if re.search(pt_e14_b23, l) is not None]

cell_meta.loc[
    cell_meta.index[row_index], ['brainregion']] = 'ERC'
cell_meta.loc[
    cell_meta.index[row_index], ['sex']] = 'Male'
cell_meta.loc[
    cell_meta.index[row_index], ['rep']] = 'MaleA'
cell_meta.loc[
    cell_meta.index[row_index], ['modality']] = 'H3K9me3'

pt_e14_b24 = "Exp14,\w\w\w:\w\w:\w\w:24"
row_index = [i-1
    for i,l in enumerate(metalines)
    if re.search(pt_e14_b24, l) is not None]

cell_meta.loc[
    cell_meta.index[row_index], ['brainregion']] = 'AMY'
cell_meta.loc[
    cell_meta.index[row_index], ['sex']] = 'Male'
cell_meta.loc[
    cell_meta.index[row_index], ['rep']] = 'MaleA'
cell_meta.loc[
    cell_meta.index[row_index], ['modality']] = 'H3K9me3'

# turn out that they are right.
# see how we generate the file in 01.clustering/src/main/R/01.read10X.R

