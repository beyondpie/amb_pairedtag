"""
Perform downsampling on Paired-Tag data and Allen Data for integration.
"""
import os
from typing import Dict, List
from subprocess import Popen
import pandas as pd
import anndata as ad
import snapatac2 as sa2
import tmpPypkg.globalvar as gv
from tmpPypkg.singlecell import downsample_ann

# * meta
projd = gv.pt_projd
cm_fnm = os.path.join(projd, "meta", "pt.barcode.meta.L5.csv")
sa2d = os.path.join(projd, "data", "pairedtag_ann")
sa2ann_fnm = os.path.join(sa2d, "pt.RNAseq.sa2ann4L1.h5ad")
k8_fnm = os.path.join(projd, "meta", "AIT21_k8_markers.txt")
# conda_dir = "/tscc/nfs/home/szu/miniforge3/"
conda_dir = "//home/szu/miniforge3/"
conda_bin = os.path.join(conda_dir, "bin/conda")
r_env = "r"
rbin = os.path.join(conda_dir, "envs", r_env, "bin", "Rscript")
sa2_env = "sa2"

tool_dir = os.path.join(projd, "extratools", "singlecell")
sa2ann2seurat_script = os.path.join(tool_dir,
                                    "SnapATACAnnData2Seurat.R")
ptd = os.path.join(gv.pt_projd, "data/pairedtag_ann")

# * add Zhaoning's annotation on L1-level.
L1toclass: Dict[int, str] = {
    0: "Neu",
    1: "Neu",
    2: "NN",
    3: "Neu",
    4: "NN",
    5: "Neu",
    6: "Neu",
    7: "NN",
    8: "Neu",
    9: "Neu",
    10: "Neu",
    11: "NN",
    12: "Neu",
    13: "NN",
    14: "NN",
    15: "NN",
    16: "Neu"
}
L1_NN: List[int] = []
L1_Neu: List[int] = []
for k, v in L1toclass.items():
    if v == "NN":
        L1_NN.append(k)
    else:
        L1_Neu.append(k)
cellmeta: pd.DataFrame = pd.read_csv(cm_fnm, header=0, sep=',')
cellmeta_NN = cellmeta[cellmeta.L1.isin(L1_NN)]
cellmeta_Neu = cellmeta[cellmeta.L1.isin(L1_Neu)]

# * load SnapATAC2 h5ad
sa2ann = sa2.read(sa2ann_fnm, backed='r')
sa2ann_NN = sa2ann.subset(
    obs_indices=cellmeta_NN.barcode.tolist(),
    out=os.path.join(sa2d, "pt.RNAseq.sa2ann.NN.h5ad")
)

fnm = os.path.join(sa2d, "pt.RNAseq.sa2ann.Neu.h5ad")
sa2ann_Neu = sa2ann.subset(
    obs_indices=cellmeta_Neu.barcode.tolist(),
    out=os.path.join(sa2d, "pt.RNAseq.sa2ann.Neu.h5ad")
)

# * further save k8 features
with open(k8_fnm, 'r') as f:
    k8s: List[str] = [g.strip() for g in f.readlines()]
k8_sa2ann: List[bool] = [
    True if g in k8s else False for g in sa2ann.var_names
]

sa2ann_NN_k8 = sa2ann_NN.subset(
    var_indices=k8_sa2ann,
    out=os.path.join(sa2d, "pt.RNAseq.sa2ann.NN.k8.h5ad")
)

sa2ann_Neu_k8 = sa2ann_Neu.subset(
    var_indices=k8_sa2ann,
    out=os.path.join(sa2d, "pt.RNAseq.sa2ann.Neu.k8.h5ad")
)

sa2ann_k8 = sa2ann.subset(
    var_indices=k8_sa2ann,
    out=os.path.join(sa2d, "pt.RNAseq.sa2ann.k8.h5ad")
)

# * sa2ann to Seurat
def get_ann2seurat_exp(sa2annfnm: str,
                       outfnm: str,
                       nds: int = 50,
                       dscol: str = "L1_2_3_4_5") -> str:
    r = [rbin, sa2ann2seurat_script, f"-f {sa2annfnm}",
         f"-m {cm_fnm}",
         f"-o {outfnm}", f"--conda {conda_bin}",
         f"--condaenv {sa2_env}", "-d", f"--dscol {dscol}",
         f"--nds {nds}"]
    return ' '.join(r)


outdir = os.path.join(projd, "data", "pairedtag_seurat")

nds = 50
cmd = get_ann2seurat_exp(
    sa2annfnm=os.path.join(sa2d, "pt.RNAseq.sa2ann.Neu.k8.h5ad"),
    outfnm=os.path.join(outdir, f"ptRNA.neu.k8.L5ds{nds}.rds"),
    nds=nds,
    dscol="L1_2_3_4_5"
)

p = Popen(cmd, shell=True)
p.wait()

# * check ds neu
a = sa2.read(os.path.join(sa2d, "pt.RNAseq.sa2ann.Neu.k8.h5ad"),
             backed='r')
cm = pd.read_csv(cm_fnm, header=0, sep=",")

barcodes_nn = cm.barcode[cm.L1.isin(L1_NN)]
barcodes_neu = cm.barcode[cm.L1.isin(L1_Neu)]
pd.Series(a.obs_names).isin(barcodes_nn).sum()
pd.Series(a.obs_names).isin(barcodes_neu).sum()

cellmeta_Neu.barcode.isin(barcodes_nn).sum()
cellmeta_Neu.barcode.isin(barcodes_neu).sum()

# * downsample nn and to Seurat
sa2ann_NN_k8 = sa2.read(
    os.path.join(sa2d, "pt.RNAseq.sa2ann.NN.k8.h5ad"), backed="r"
)
# check all the barcodes in cellmeta_NN
barcodes_NN = cellmeta_NN.barcode
pd.Series(sa2ann_NN_k8.obs_names).isin(barcodes_NN).sum()

nds = 30
cmd = get_ann2seurat_exp(
    sa2annfnm=os.path.join(sa2d, "pt.RNAseq.sa2ann.NN.k8.h5ad"),
    outfnm=os.path.join(outdir, f"ptRNA.nn.k8.L5ds{nds}.rds"),
    nds=nds,
    dscol="L1_2_3_4_5"
)

p = Popen(cmd, shell=True)
p.wait()

# * partition into region-specific neuronal cells
sa2_neu_region_dir = os.path.join(
    sa2d, "neu_region")
os.makedirs(sa2_neu_region_dir, exist_ok=True)
cm_fnm = os.path.join(projd, "meta",
                      "pairedtag.cell.meta.all.with.init.tf.csv")
cellmeta: pd.DataFrame = pd.read_csv(cm_fnm, header=0, sep=',')
sa2_neu: sa2.AnnData = sa2.read(
    filename=os.path.join(sa2d, "pt.RNAseq.sa2ann.Neu.k8.h5ad"),
    backed='r'
)

allen2pt_region: Dict[str, List[str]] = {
    "AMY": ["AMY"],
    "CPU": ["CPU"],
    "ERC": ["ERC"],
    "HIP": ["HCa", "HCp"],
    "HYP": ["HYP"],
    "NAC": ["NAC"],
    "PFC": ["PFC"],
    "VTA": ["VTA_SnR"]
}


def partitionByRegion(r: str) -> sa2.AnnData:
    outfnm = os.path.join(sa2_neu_region_dir,
                          f"pt.neu.{r}.all.sa2.h5ad")
    barcodes = cellmeta[cellmeta.brainregion.isin(
        allen2pt_region[r])].barcode
    barcodes = barcodes[barcodes.isin(sa2_neu.obs_names)]
    print(f"region {r}: {len(barcodes)}.")
    print(f"outfile to: {outfnm}")
    r = sa2_neu.subset(obs_indices=barcodes, out=outfnm)
    return r


regions = list(allen2pt_region.keys())
regions.remove("AMY")
for r in regions:
    partitionByRegion(r=r)

# * transform the region-specific anndata to seurat object
sa2_neu_region_dir = os.path.join(
    sa2d, "neu_region")
cm_fnm = os.path.join(projd, "meta",
                      "pairedtag.cell.meta.all.with.init.tf.csv")
sa2_neu_region_seu_dir = os.path.join(
    projd, "data", "pairedtag_seurat", "neu_seu_region")
os.makedirs(sa2_neu_region_seu_dir, exist_ok=True)
def get_neuann2seu_exp(r: str) -> str:
    sa2annfnm = os.path.join(sa2_neu_region_dir,
                             f"pt.neu.{r}.all.sa2.h5ad")
    outseufnm = os.path.join(sa2_neu_region_seu_dir,
                             f"pt.neu.{r}.all.seu.rds")
    r = [rbin, sa2ann2seurat_script, f"-f {sa2annfnm}",
         f"-o {outseufnm}", f"--conda {conda_bin}", f"-m {cm_fnm}",
         f"--condaenv {sa2_env}"]
    return ' '.join(r)


mrs = list(allen2pt_region.keys())
for r in mrs:
    cmd = get_neuann2seu_exp(r=r)
    p = Popen(cmd, shell=True)
    p.wait()

# * regenearate down sampling based on our latest annotation
# 2024-10-05
cellMeta = pd.read_csv(gv.ptcellmetafnm, sep=",", header=0)
cellMeta = cellMeta.loc[cellMeta.annotQuality == "Good", ]
cellMeta.set_index('barcode', drop=False, inplace=True)

class_cells = cellMeta["annot.sc"].apply(lambda x: x.split(" ")[-1])

nn_cells = cellMeta.loc[class_cells == "NN", ].barcode
neu_cells = cellMeta.barcode[~cellMeta.barcode.isin(nn_cells)]
ann: ad.AnnData = sa2.read(
    os.path.join(ptd, "pt.RNAseq.sa2ann.k8.h5ad"), backed=None)
ann.obs.insert(0, 'barcode', ann.obs_names)
ann = ann[cellMeta.index, ].copy()
ann.obs.insert(2, 'subclass', cellMeta.loc[ann.obs_names, "annot.sc"])

ann_nn = ann[nn_cells, ].copy()
ann_neu = ann[neu_cells, ].copy()
ann_nn.write_h5ad(os.path.join(ptd, "pt.RNA.ann.nn.noimn.k8.h5ad"))
ann_neu.write_h5ad(os.path.join(ptd, "pt.RNA.ann.neu.withimn.k8.h5ad"))

# perform downsample
ann_nnds = downsample_ann(ann_nn, on="subclass", nds=10000)
ann_neuds = downsample_ann(ann_neu, on="subclass", nds=1000)
ann_nnds.write_h5ad(
    os.path.join(ptd, "pt.RNA.ann.nn.noimn.k8.ds.h5ad"))
ann_neuds.write_h5ad(
    os.path.join(ptd, "pt.RNA.ann.neu.withimn.k8.ds.h5ad"))

# 2024-11-12
# 1.8 million cells
ann_neu = ad.read_h5ad(os.path.join(ptd, "pt.RNA.ann.neu.withimn.k8.h5ad"),
                       backed='r')

# 2024-11-18
ann_nnds = ad.read_h5ad(os.path.join(ptd, "pt.RNA.ann.nn.noimn.k8.ds.h5ad"),
                        backed=None)
ann_neuds = ad.read_h5ad(os.path.join(ptd, "pt.RNA.ann.neu.withimn.k8.ds.h5ad"),
                         backed=None)
ann_ds = ad.concat([ann_nnds, ann_neuds], merge="same")
ann_ds.write_h5ad(os.path.join(ptd, "pt.RNA.ann.k8.ds.raw.h5ad"))

