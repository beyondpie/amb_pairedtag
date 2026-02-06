"""
- Calculate the similarities between subclasses based on
  four different histone modifications.

- H3K9me3: find the consistent repressive regions.
  - used for scatter plot of active signals vs repressive signals
  - used for filter out the regions we are not interested for
    correlation and other analysis
"""
import os
from typing import Dict
import numpy as np
import pandas as pd

import snapatac2 as sa2
import anndata as ad

import matplotlib.pyplot as plt
import seaborn as sns

import tmpPypkg.globalvar as gv
from tmpPypkg.globalvar import reorder_scs
from tmpPypkg.globalvar import transform_allenlabel
from tmpPypkg.bam import call_plotcor

# * meta
avgPeakSize: Dict[str, int] = {
    "H3K27ac": 500,
    "H3K4me1": 1200,
    "H3K27me3": 1000,
    "H3K9me3": 700
}

# To identify the consitutive repressive region, 
# H3K9me3 size is set as 10kb
# (smaller region could be better but let's use this firstly)

annd = os.path.join(gv.pt_projd, "data", "pairedtag_ann")
h3k4me1_ann_fnm = os.path.join(annd,
                               "2024701.brain.snapatac2.h3k4me1.h5ad")
h3k27ac_ann_fnm = os.path.join(annd,
                               "2024701.brain.snapatac2.h3k27ac.h5ad")
h3k27me3_ann_fnm = os.path.join(annd,
                                "20240821.ambPT.snapatac2.h3k27me3.10kbin.h5ad")
h3k9me3_ann_fnm = os.path.join(annd,
                               "20240821.ambPT.snapatac2.h3k9me3.10kbin.h5ad")

# * organize AnnData with cellMeta
hcellMeta: pd.DataFrame = gv.getPairedTagCellMeta()


# * load AnnData
# obs is pandas's DataFrame
h3k9me3_ann = sa2.read(filename=h3k9me3_ann_fnm, backed=None)
h3k4me1_ann = sa2.read(filename=h3k4me1_ann_fnm, backed=None)


# * get pseudobulk-level expressions
h3k9me3_barcodes = hcellMeta.barcode
h3k9me3_ann = h3k9me3_ann[
    h3k9me3_ann.obs_names.isin(h3k9me3_barcodes)]
h3k9me3_sc = sa2.tl.aggregate_X(
    adata=h3k9me3_ann, groupby="annot.sc", normalize="RPM")
h3k9me3_sc.write_h5ad(
    filename=os.path.join(annd, "20240821.h3k9me3.ann.10kb.pseudosc.RPM.h5ad"))


# * select variable features for H3K9me3

# * H3K27me3 not on stable H3K9me3 regions

# * active enhancer: H3K27ac + H3K4me1

# * poised enhancer: -H3K27ac + H3K4me1
