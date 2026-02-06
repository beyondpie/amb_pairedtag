import sys
import os
import snapatac2 as sa2
import numpy as np
import pandas as pd
import polars as pl
import anndata as ad
import scanpy as sc
import importlib
importlib.reload(sc)


# * meta
projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
annd = os.path.join(projd, "data", "pairedtag_ann")
workd = os.path.join(projd, "17.repressiveMarks")
outd = os.path.join(workd, "out", "allGenomicRange")

# * functions
# * main
annK27me3 = sa2.read(
    filename = os.path.join(annd,
                            "2024701.brain.snapatac2.h3k27me3.50kbin.h5ad"),
    backed = None
)

annK9me3 = sa2.read(
    filename = os.path.join(annd,
                            "2024701.brain.snapatac2.h3k9me3.50kbin.h5ad"),
    backed = None
)

# output raw CPM on selected top features

## somehow two exp columns in obs
annK27me3.obs.drop("exp", axis =1, inplace = True)
cell_barcode = annK27me3.obs.index[annK27me3.obs["annotQuality"] == "Good"]
feas  = annK27me3.var.index[annK27me3.var["selected"]]
sK27me3 = annK27me3[cell_barcode.to_list(), feas.to_list()]
aggK27me3 = sc.get.aggregate(sK27me3, by = "annot.sc", func = "sum")
rpmK27me3 = pd.DataFrame(aggK27me3.layers["sum"],
                         index = aggK27me3.obs["annot.sc"],
                         columns = aggK27me3.var_names)
rpmK27me3.to_csv(os.path.join(outd, "aggK27me3.pd.csv"), sep = ",",
                 header = True, index = True)


## somehow two exp columns in obs
annK9me3.obs.drop("exp", axis =1, inplace = True)
cell_barcode = annK9me3.obs.index[annK9me3.obs["annotQuality"] == "Good"]
feas  = annK9me3.var.index[annK9me3.var["selected"]]
sK9me3 = annK9me3[cell_barcode.to_list(), feas.to_list()]
aggK9me3 = sc.get.aggregate(sK9me3, by = "annot.sc", func = "sum")
rpmK9me3 = pd.DataFrame(aggK9me3.layers["sum"],
                         index = aggK9me3.obs["annot.sc"],
                         columns = aggK9me3.var_names)
rpmK9me3.to_csv(os.path.join(outd, "aggK9me3.pd.csv"), sep = ",",
                 header = True, index = True)

## on class level
aggClassK27me3 = sc.get.aggregate(sK27me3, by = "annot.c", func = "sum")
rpmClassK27me3 = pd.DataFrame(
    aggClassK27me3.layers["sum"],
    index = aggClassK27me3.obs["annot.c"],
    columns = aggClassK27me3.var_names
    )
rpmClassK27me3.to_csv(os.path.join(outd, "aggClassK27me3.pd.csv"),
                      sep = ",", header = True, index = True)


aggClassK9me3 = sc.get.aggregate(sK9me3, by = "annot.c", func = "sum")
rpmClassK9me3 = pd.DataFrame(
    aggClassK9me3.layers["sum"],
    index = aggClassK9me3.obs["annot.c"],
    columns = aggClassK9me3.var_names
    )
rpmClassK9me3.to_csv(os.path.join(outd, "aggClassK9me3.pd.csv"),
                      sep = ",", header = True, index = True)
