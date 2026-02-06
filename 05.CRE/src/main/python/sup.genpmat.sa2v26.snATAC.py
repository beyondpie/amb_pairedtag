from typing import List
import pandas as pd
import os
import sys
import snapatac2 as sa2
import re

peakbedfnm = sys.argv[1]
tag = sys.argv[2]

if not os.path.exists(peakbedfnm):
    raise RuntimeError(f"{peakbedfnm} does not found.")

# * meta
projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
outd = os.path.join(projd, "data", "snATAC")
outfnm = os.path.join(outd, f"all.sa2v26.pmat.{tag}.nofrag.h5ad")

# * functions
def transform_allenlabel(a: str) -> str:
    return re.sub(r" +|/|-", "_", a)


# * main
print("Read anndataset.")
all_adset = sa2.read_dataset(
    os.path.join(projd, "data", "snATAC",
                 "all.sa2v26.adset.h5ad"), mode="r"
)

print(f"Start to add pmat for {peakbedfnm} ...")
annpmat = sa2.pp.make_peak_matrix(
    adata=all_adset,
    peak_file=peakbedfnm,
    inplace=False,
    file=outfnm,
    chunk_size=5000
)

all_adset.close()
print(f"Generate the pmat anndata: {outfnm}.")

# * add cell meta to the pmats
print("Add subclass, cluster and ptgroup info.")

snATAC_cellmeta = pd.read_csv(
   os.path.join("/projects/ps-renlab2/szu/projects/amb_pairedtag",
                "data", "snATAC", "snATAC_cellmeta.v9.8.240830.tsv"),
    sep="\t", header=0
)
snATAC_cellmeta.set_index("barcode2", drop=False, inplace=True)

def add_col(ann: sa2.AnnData) -> None:
    ann_barcode: List[str] = ann.obs_names
    cm = snATAC_cellmeta.loc[ann_barcode]
    ann.obs["subclass"] = cm[
        "subclass_id_label_v3"].apply(transform_allenlabel).to_list()
    ann.obs["cluster"] = cm["L4"].to_list()
    ann.obs["ptgroup"] = cm["ptgroup"].to_list()


add_col(annpmat)

annpmat.close()
print("Done.")
