"""
1. get sublibrary_id to bam files
2. get cell barcode to out bam dir and fnm
"""

import sys
import os
from typing import List, Dict, Tuple
from collections import OrderedDict

import pandas as pd
import pyprojroot

proj_dir = pyprojroot.here()
meta_sublib = pd.read_csv(
    os.path.join(proj_dir, "meta", "amb_pairedtag_sublibrary.csv"),
    sep = ",",
    header = 0,
    names = None
)


def get_mapping_dir(path: str) -> str:
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    dir1 = os.path.join(path, "03.mapping")
    if os.path.exists(dir1):
        return dir1
    dir2 = os.path.join(path, "03.mm10_mapping")
    if os.path.exists(dir2):
        return dir2
    raise FileNotFoundError(f"No mapping dir: {path}")

map_dir = meta_sublib.apply(
    lambda x: get_mapping_dir(x['path']), axis = 1).to_list()
        
sub_DNAs = meta_sublib['DNA'].to_list()

mapping_fnms = [os.path.join(i, f"{j}_mm10_sorted_rmdup.bam") for i, j in
                zip(map_dir, sub_DNAs)]

# check exist
sum([os.path.exists(i) for i in mapping_fnms])  == len(map_dir)

meta_sublib['DNA_bam'] = mapping_fnms

meta_sublib.to_csv(
    os.path.join(proj_dir, "meta", "amb_pairedtag_sublibrary_v2.csv"),
    sep = ",",
    header = True,
    index = False
)

# * add RNA bam files to amb pairedtag sublibrary
# 20240321
meta_sublib = pd.read_csv(
    os.path.join(proj_dir, "meta", "amb_pairedtag_sublibrary_v2.csv"),
    sep = ",",
    header = 0,
    names = None
)

RNA_bam_fnms = meta_sublib.apply(
    lambda x: x['DNA_bam'].replace(x['DNA'], x['RNA']), axis = 1)
# check file exist
sum([os.path.exists(i) for i in RNA_bam_fnms]) == meta_sublib.shape[0]

meta_sublib['RNA_bam'] = RNA_bam_fnms
# here we remove the version label
meta_sublib.to_csv(
    os.path.join(proj_dir, "meta", "amb_pairedtag_sublibrary.csv"),
    sep = ",",
    header = True,
    index = False
)

# * get cell barcode to out bam dir
cell_meta_fnm = os.path.join(proj_dir,
                                "meta", "pairedtag.meta.v1.csv")
cell_meta = pd.read_csv(cell_meta_fnm,
                        sep = ",", header = 0,
                        names = None)
barcode2out_bam_dir: List[Tuple[str, str]] = [
    (row['barcode'] ,
      "{r}-{h}-{s}".format(
          r = row['brainregion'],
          h = row['modality'],
          s = row['rep']))
     for _, row in cell_meta.iterrows()]
with open(os.path.join(proj_dir, "meta", "barcode2bamdir.csv"), 'w') as f:
    f.write("barcode,bamdir\n")
    f.writelines([f"{i[0]},{i[1]}\n" for i in barcode2out_bam_dir])



