import os
from typing import Dict, List
import pandas as pd
import scanpy
import snapatac2 as sa2
import re

# * functions
def transform_allenlabel(a: str) -> str:
    return re.sub(r" +|/|-", "_", a)


# * meta
projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"

snATAC_cellmeta_fnm = os.path.join(
    "/projects/ps-renlab2/szu/projects/CEMBA2", "meta/wmb.cellmeta.v9.7.tsv"
)

snATAC_fragd = os.path.join(
    "/projects/ps-renlab2/szu/projects/CEMBA2",
    "00.data.preprocess",
    "snapatac2_pp_out/raw_fragment_files",
)

snATAC_bamd = os.path.join("/projects/ps-renlab/szu/shared_data", "wmb_snATAC/raw_bams")

mm10_chrom_fnm = os.path.join(projd, "meta", "mm10.chrom.sizes.lite")

snATAC_samples = os.path.join(
    projd,
    "data/snATAC",
    "snapatac2.6_samples",
)

blacklist_fnm = os.path.join(
    "/projects/ps-renlab2/szu/projects/amb_pairedtag", "meta", "mm10-blacklist.v2.bed"
)

snATAC_sa2d = os.path.join(
    "/projects/ps-renlab2/szu/projects/amb_pairedtag", "data", "snATAC"
)

tmpdir = "/projects/ps-renlab2/szu/tmpdir"

# * main
snATAC_cellmeta = pd.read_csv(snATAC_cellmeta_fnm, sep="\t", header=0)

fragfnms = [
    f for f in os.listdir(snATAC_fragd) if os.path.isfile(os.path.join(snATAC_fragd, f))
]

all_samples = [f.split(".")[0] for f in fragfnms]
with open(
    os.path.join(projd, "05.CRE/src/main/resource", "snATACsample.txt"), "w"
) as f:
    f.write("\n".join(all_samples))

for s in all_samples:
    barcodes = snATAC_cellmeta[snATAC_cellmeta["sample"] == s].barcode.to_list()
    outf = os.path.join(
        projd, "05.CRE/src/main/resource", "snATAC_sample2barcode", f"{s}.sample.txt"
    )
    with open(outf, "w") as f:
        f.write("\n".join(barcodes))
        print(f"write {len(barcodes)} to {outf}.")

# * merge all the h5ads
H3K27ac_peak_fnm = os.path.join(
    "/projects/ps-renlab2/zhw063/99.MouseBrainPairedTag/",
    "snapatac_dna_mod_analysis",
    "20240822_marker_peaks_H3K27ac",
    "subclass_heatmap_peaks_metadata.csv",
)

name2peak = pd.read_csv(H3K27ac_peak_fnm, sep=",", header=0)
peaks = name2peak["Genomic Location"].apply(
    lambda x: "\t".join([x.replace(":", "\t").replace("-", "\t"), x])
)
peak_bed_fnm = os.path.join(
    "/projects/ps-renlab2/szu/projects",
    "amb_pairedtag",
    "05.CRE/src/main/resource",
    "subclas_heatmap_peaks_metadata.csv",
)

peaks.to_csv(peak_bed_fnm, header=False, index=False)

# get anndataset
with open(
    os.path.join(projd, "05.CRE/src/main/resource", "snATACsample.txt"), "r"
) as f:
    all_samples = [s.strip() for s in f.readlines()]
s2anns = [
    (s, sa2.read(filename=os.path.join(snATAC_samples, f"{s}.sa2v26.h5ad"), backed="r"))
    for s in all_samples
]


all_adset = sa2.AnnDataSet(
    adatas=s2anns,
    filename=os.path.join(snATAC_sa2d, "all.sa2v26.adset.h5ad"),
    add_key="sample",
)

# update barcode name
old_obsnames = all_adset.obs_names
sample_info = all_adset.obs["sample"]
new_obsnames = [f"{s}.{b}" for s, b in zip(sample_info, old_obsnames)]
all_adset.obs_names = new_obsnames
all_adset.close()

# add H3K27ac-based peaks
all_adset = sa2.read_dataset(
    filename=os.path.join(snATAC_sa2d, "all.sa2v26.adset.h5ad"), mode="r"
)

print("start to add pmat.")
all_adset_pmat = sa2.pp.make_peak_matrix(
    adata=all_adset,
    peak_file=peak_bed_fnm,
    inplace=False,
    file=os.path.join(snATAC_sa2d, "all.sa2v26.ann.pmat.nofrag.h5ad"),
)
all_adset.close()
all_adset_pmat.close()
print("add pmat done.")
all_ann = all_adset.to_adata(file=os.path.join(snATAC_sa2d, "all.sa2v26.ann.h5ad"))

print("load pmat-contained ann dataset...")
all_only_pmat = sa2.read(
    filename=os.path.join(snATAC_sa2d, "all.sa2v26.ann.pmat.nofrag.h5ad"),
    backed="r"
)

# * add metadata to anndata with pmat

# snATAC raw cell meta
snATAC_cellmeta = pd.read_csv(snATAC_cellmeta_fnm, sep="\t", header=0)
snATAC_cellmeta.set_index(keys="barcode2", drop=False, inplace=True)

# nn supertype
L4annot = pd.read_csv(
    os.path.join(snATAC_sa2d, "snATACTransferLabel",
                 "snATAC_nnL4tosp_240814.csv"),
    sep=",", header=0
)
L4annot.set_index("L4", drop=False, inplace=True)

# add the sp to snATAC_cellmeta
ptgroup = snATAC_cellmeta[
    "subclass_id_label_v3"].copy().apply(transform_allenlabel)

idx_nn = snATAC_cellmeta["L4"].isin(L4annot.index)
supertype = L4annot.loc[
    snATAC_cellmeta.loc[idx_nn]["L4"]]["sp"].apply(transform_allenlabel)
ptgroup[idx_nn] = supertype
snATAC_cellmeta["ptgroup"] = ptgroup
snATAC_cellmeta.to_csv(
    os.path.join(snATAC_sa2d, "snATAC_cellmeta.v9.8.240830.tsv"),
    sep="\t",
    header=True,
    index=False
)


annpmat = sa2.read(os.path.join(snATAC_sa2d, "all.sa2v26.ann.pmat.nofrag.h5ad"),
                   backed="r+")

ann_barcodes: List[str] = annpmat.obs_names
ann_old_obs = annpmat.obs


subclass = snATAC_cellmeta.loc[
    ann_barcodes]["subclass_id_label_v3"].apply(transform_allenlabel).to_list()
cluster = snATAC_cellmeta.loc[
    ann_barcodes]["L4"].to_list()
ptgroup = snATAC_cellmeta.loc[ann_barcodes]["ptgroup"].to_list()

annpmat.obs['subclass'] = subclass
annpmat.obs['cluster'] = cluster
annpmat.obs['ptgroup'] = ptgroup
annpmat.close()


annpmat = sa2.read(os.path.join(snATAC_sa2d, "all.sa2v26.ann.pmat.nofrag.h5ad"),
                   backed="r")
