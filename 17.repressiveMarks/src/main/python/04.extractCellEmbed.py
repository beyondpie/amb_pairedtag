import numpy as np
import pandas as pd
import anndata as ad
import snapatac2 as sa2
import os

# * meta
projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
annd = os.path.join(projd, "data", "pairedtag_ann")
workd = os.path.join(projd, "17.repressiveMarks")
outd = os.path.join(workd, "out", "allGenomicRange")

# * functions
def loadsa2annToAnnMem(mod: str) -> ad.AnnData:
    """
    Returns view of original AnnData.
    Use .copy() if needing independent AnnData.
    """
    r: ad.AnnData = sa2.read(os.path.join(
        annd, f"2024701.brain.snapatac2.{mod}.50kbin.h5ad"),
             backed = None)
    r.obs.drop("exp", axis = 1, inplace = True)
    cells = r.obs.index[r.obs["annotQuality"] == "Good"]
    return r[cells, :]

def getCellEmbed(x: ad.AnnData) -> pd.DataFrame:
    r1 = x.obs[["annot.sc", "annot.c"]]
    r2 = pd.DataFrame(
        data = x.obsm["X_spectral"],
        index = x.obs_names,
        columns = [ f"S{i}"for i in range(x.obsm["X_spectral"].shape[1])]
    )
    r3 = pd.DataFrame(
        data = x.obsm["X_umap"],
        index = x.obs_names,
        columns = [f"U{i}" for i in range(2)]
    )
    # combine two data.frame
    r = pd.concat([r1, r2, r3], axis = 1)
    return r

def saveCellEmbed(x: pd.DataFrame, mod: str) -> None:
    outfnm = os.path.join(outd, f"pt.singlecell.emebd.{mod}.csv")
    x.to_csv(outfnm, sep = ",", header = True, index = True)
    return None
    
# * main
annK27me3 = loadsa2annToAnnMem(mod = "h3k27me3")
rK27me3 = getCellEmbed(annK27me3)
saveCellEmbed(rK27me3, mod = "H3K27me3")

annK9me3 = loadsa2annToAnnMem(mod = "h3k9me3")
rK9me3 = getCellEmbed(annK9me3)
saveCellEmbed(rK9me3, mod = "H3K9me3")


