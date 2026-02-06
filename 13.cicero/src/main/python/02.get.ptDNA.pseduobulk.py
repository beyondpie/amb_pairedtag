import os
import pandas as pd
import anndata as ad
import snapatac2 as sa2
import tmpPypkg.globalvar as gv
from adpbulk import ADPBulk

# * meta
projd = gv.pt_projd

# * H3K27ac
ann: ad.AnnData = sa2.read(
    os.path.join(projd, "data", "pairedtag_ann",
                 "ann.H3K27ac.pmat.h5ad"),
    backed=None)
adpb = ADPBulk(ann, "annot.sc")
pseudobulk_matrix = adpb.fit_transform()
pseudobulk_matrix.to_hdf(
    os.path.join(projd, "data", "pairedtag_ann",
                 "pt.H3K27ac.pseudobulk.count.sc.bestSPM.h5"),
    key="H3K27ac", mode="w", index=True
)
