import os
import numpy as np
import anndata as ad
import snapatac2 as sa2
import tmpPypkg.globalvar as gv


projd = gv.pt_projd

ann = ad.read_h5ad(
    os.path.join(projd, "12.DE", "out",
                 "ann.H3K27ac.pmat.h5ad"), backed=None)
sa2ann = sa2.read(
    os.path.join(projd, "data", "pairedtag_ann", "2024701.brain.snapatac2.h3k27ac.h5ad"),
    backed="r"
)

sa2ann = sa2.read(os.path.join(projd, "12.DE", "out", "ann.H3K27ac.pmat.h5ad"), backed="r")


a = sa2ann.subset(obs_indices=sa2ann.obs_names[0:100], inplace=False,
                  out=None,
                  backend="hdf5")

fea = sa2.pp.select_features(
    adata=a,
    n_features=50000,
    inplace=False
)


