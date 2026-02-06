import os
import pandas as pd
import anndata as ad
import snapatac2 as sa2
import tmpPypkg.globalvar as gv
from adpbulk import ADPBulk

# * meta
projd = gv.pt_projd
ann: ad.AnnData = sa2.read(
    os.path.join(projd, "data", "pairedtag_ann",
                 "pt.RNAseq.sa2ann.filterbymito.h5ad"),
    backed=None
)

cellmeta: pd.DataFrame = pd.read_csv(
    os.path.join(projd, "meta", "pairedtag.cell.meta.all.240626.csv"),
    header=0, sep=",")
cellmetaQC: pd.DataFrame = cellmeta[cellmeta.annotQuality == "Good"]
cellmetaQC.set_index("barcode", inplace=True)

annQC = ann[ann.obs_names.isin(cellmetaQC.barcode)].copy()
annQC_obs = cellmetaQC.loc[annQC.obs_names]
annQC.obs = annQC_obs

annQC.write_h5ad(filename=os.path.join(
    projd, "data", "pairedtag_ann",
                 "pt.RNAseq.ann.goodQuality.fullmeta.allgenesymbol.h5ad"))

# * get pseudo-bulk expression
adpb = ADPBulk(annQC, "annot.sc")
pseudobulk_matrix = adpb.fit_transform()
sample_meta = adpb.get_meta()
pseudobulk_matrix.to_csv(
    os.path.join(projd, "data", "pairedtag_ann",
                 "pt.RNAseq.pseudobulk.count.sc.allgenesymbol.pdDataFrame.csv"),
    sep=",",
    header=True, index=True
)
