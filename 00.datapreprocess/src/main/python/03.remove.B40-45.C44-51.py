# Cluster 11 cells are mainly coming from early
# library preps of Exp2 (B40-B45) and Exp3 (C44-C51), where we had a suboptimal batch of
# library prep primers. There does not seem to be a brain-region bias where these cells
# come from. I did stand-alone analysis on these sequencing samples, which looks fine.
# However, they carry some low-complexity sequences due to that primer issue, which can
# be mapped to transcripts like Gm10800 and Gm10801, which are T-rich low complexity
# sequences.
# so we now remove the B40-45 and C44-C51.

import os
import sys
from typing import List
import pandas as pd
import pyprojroot

proj_root:str = pyprojroot.here()
cell_meta :pd.DataFrame = pd.read_csv(
    os.path.join(proj_root, "meta", "pairedtag.meta.v1.csv"),
    sep = ",",
    header = 0
)
cell_meta.set_index("barcode", drop = False, inplace = True)
all_barcode = cell_meta.barcode
cell_meta.drop(columns="barcode", inplace = True)
cell_meta.insert(0, "barcode", all_barcode)

ncol: int = cell_meta.shape[1]

cell_dlt: pd.DataFrame = pd.read_csv(
    os.path.join(proj_root, "meta", "sa2_scrublet_8k_all.csv"),
    sep = ",",
    header = 0
)
cell_dlt.set_index("barcode", drop = False, inplace = True)

cell_meta.insert(ncol, "dltprob", cell_dlt.loc[cell_meta.barcode, "dlt_prob"])
cell_meta.insert(ncol+1, "dltscore", cell_dlt.loc[cell_meta.barcode, "dlt_score"])

barcode2dlt_afdlt: pd.DataFrame = pd.read_csv(
    os.path.join(proj_root, "meta", "barcode2dlt.filtered.L1.k8.dltscore.csv"),
    sep = ",",
    header = 0
)
barcode_afdlt = set(barcode2dlt_afdlt.barcode)
isdlt: List[str] = [
    "no" if b in barcode_afdlt else "yes" for b in cell_meta.barcode]
cell_meta.insert(ncol+2, "isdlt", isdlt)

libs = cell_meta.barcode.apply(lambda x : x.split(":")[0] )
bias_libs: List[str] = [
    f"B{i}" for i in range(40, 46)] + [f"C{i}" for i in range(44,52)]
is_biaslib: List[str] = [
    "no" if l not in bias_libs else "yes" for l in libs
]
cell_meta.insert(ncol+3, "isbiaslib", is_biaslib)

isLQ = (cell_meta.isdlt == "yes") | (cell_meta.isbiaslib == "yes")

cell_meta.insert(ncol + 4, "isLQ", isLQ)

cell_meta.to_csv(
    os.path.join(proj_root, "meta", "pt.barcode.meta.with.qc.v2.csv"),
    sep = ",",
    header = True,
    index = False
)





