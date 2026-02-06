import os
import re
from typing import List, Dict
import functools
from subprocess import Popen
import pandas as pd
import anndata as ad

import tmpPypkg.globalvar as gv

# * configs
allen_dir = os.path.join("/projects/ps-renlab2/szu/projects", "shared_data/allen")
allend = os.path.join(gv.pt_projd, "data/allen")
allen_10xv3_path = os.path.join(allen_dir, "AIT21_10Xv3.h5ad")
allen_10xv2_path = os.path.join(allen_dir, "AIT21_10Xv2.h5ad")
allen_snRNA_path = os.path.join(allen_dir, "AIT21_multiome.h5ad")

proj_root = gv.pt_projd
out_dir = os.path.join(proj_root, "data/allen")

# * save allen cell-level meta
allen_10xv3 = ad.read_h5ad(filename=allen_10xv3_path, backed="r")
ameta = allen_10xv3.obs
ameta.insert(0, "barcode", allen_10xv3.obs_names)
ameta.to_csv(
    os.path.join(out_dir, "allen.10xv3.cell.meta.csv"),
    sep=",",
    header=True,
    index=False,
)

# * split data based on region information.
# need 20 minutes to load
print(f"Loading data: {allen_10xv3_path}.")
# load to memory
# if backed = 'r', will lead later allen_10xv3 change of filename
# after write_h5ad
# and if we set original filename back, it will change our local h5ad
allen_10xv3 = ad.read_h5ad(filename=allen_10xv3_path, backed=None)
r2cnt: pd.Series = allen_10xv3.obs["roi"].value_counts()
r2cnt.to_csv(
    os.path.join(out_dir, "allen_10xv3_region2count.csv"),
    sep=",",
    header=False,
    index=True,
)
# for each region, we write it into disk
raw_regions: List[str] = r2cnt.index.to_list()
rename_regions = [re.sub(r" - ", "_", i) for i in raw_regions]
rename_regions = [re.sub(r" ", "_", i) for i in rename_regions]
raw2fnm: Dict[str, str] = {
    k: os.path.join(out_dir, f"allen_10xv3_{v}.h5ad")
    for k, v in zip(raw_regions, rename_regions)
}
for raw_region in raw_regions:
    print(f"subset region: {raw_region}...")
    tmp = allen_10xv3[allen_10xv3.obs["roi"] == raw_region]
    tmp.write_h5ad(filename=raw2fnm[raw_region])

# allen snRNA
print(f"Loading data: {allen_snRNA_path}.")
allen_snRNA = ad.read_h5ad(filename=allen_snRNA_path, backed=None)

r2cnt: pd.Series = allen_snRNA.obs["roi"].value_counts()
r2cnt.to_csv(
    os.path.join(out_dir, "allen_snRNA_region2count.csv"),
    sep=",",
    header=False,
    index=True,
)
# for each region, we write it into disk
raw_regions: List[str] = r2cnt.index.to_list()
rename_regions = [re.sub(r" - ", "_", i) for i in raw_regions]
rename_regions = [re.sub(r" ", "_", i) for i in rename_regions]
raw2fnm: Dict[str, str] = {
    k: os.path.join(out_dir, f"allen_snRNA_{v}.h5ad")
    for k, v in zip(raw_regions, rename_regions)
}
for raw_region in raw_regions:
    print(f"subset region: {raw_region}...")
    tmp = allen_snRNA[allen_snRNA.obs["roi"] == raw_region]
    tmp.write_h5ad(filename=raw2fnm[raw_region])

# * generate test h5ad for toSeurat script
# allen_dir = "/projects/ps-renlab2/szu/projects/amb_pairedtag/data/allen"
# ann = ad.read_h5ad(
#     filename=os.path.join(
#         allen_dir, "allen_snRNA_Mouse_Multiome_MY_AP-DCO-VCO-DCN-ECU-NTS.h5ad"
#     ),
#     backed=None,
# )
# target_dir = os.path.join(proj_root, "03.integration", "src/test/resource")
# ann.write_h5ad(os.path.join(target_dir, "scanpy_ann.h5ad"))
# cell_meta = ann.obs
# cell_meta.insert(0, "barcode", ann.obs_names)
# cell_meta.to_csv(
#     os.path.join(target_dir, "scanpy_ann.cellmeta.csv"),
#     sep=",",
#     header=True,
#     index=False,
# )

# * get 10xv3 allen with the regions interested.
# downsample based on subclass level
proj_root = gv.pt_projd
work_dir = os.path.join(proj_root, "03.integration")
allen_dir = os.path.join(proj_root, "data/allen")
pt2allen: pd.DataFrame = pd.read_csv(
    os.path.join(
        work_dir, "src/main/resource", "Allen_reference_data_locator_ZW_20240131.csv"
    ),
    sep=",",
    header=0,
)

allen_10xv3_fnms: List[str] = (
    pd.Series(
        functools.reduce(
            lambda a, b: a + b,
            pt2allen.iloc[:, 4]
            .apply(lambda x: x.replace(" ", "").split(";"))
            .to_list(),
        )
    )
    .unique()
    .tolist()
)

allen_10xv3_h5ads: List[ad.AnnData] = [
    ad.read_h5ad(os.path.join(allen_dir, i)) for i in allen_10xv3_fnms
]

allen_10xv3_h5ad: ad.AnnData = ad.concat(adatas=allen_10xv3_h5ads, axis=0)
# allen_10xv3_h5ad.write_h5ad(
#    os.path.join(allen_dir, "allen.10xv3.pt.regions.ann.h5ad"))
allen_meta: pd.DataFrame = pd.read_csv(
    os.path.join(work_dir, "src/main/resource", "AIT21_annotation_freeze_081523.tsv"),
    sep="\t",
    header=0,
)

allen_meta.set_index("cl", inplace=True, drop=False)
cl2subclass: Dict[int, str] = {
    row["cl"]: row["subclass_id_label"] for _, row in allen_meta.iterrows()
}
cl2supertype: Dict[int, str] = {
    row["cl"]: row["supertype_id_label"] for _, row in allen_meta.iterrows()
}
cl2ntype: Dict[int, str] = {
    row["cl"]: row["nt_type_label"] for _, row in allen_meta.iterrows()
}

suffix: pd.Series = allen_meta.subclass_label.apply(lambda x: x.split(" ")[-1])
suffix.value_counts(dropna=False)

neu_cls: List[int] = allen_meta[~suffix.str.contains("NN|LQ", regex=True)].cl.to_list()
neu_sc: pd.DataFrame = pd.DataFrame(
    [(k, cl2subclass[k]) for k in neu_cls], columns=["cl", "subclass"]
)
neu_sc.insert(column="nt_type_label", value=[cl2ntype[k] for k in neu_cls], loc=2)

# check neu_sc accuracy
# a = neu_sc.apply(
#     lambda x: f"{x.subclass}:{x.nt_type_label}", axis = 1
# ).value_counts(dropna = False)
# a.to_csv("tmp.csv")

allen_10xv3_cl_neu = allen_10xv3_h5ad.obs.cl.apply(lambda x: int(x)).isin(neu_cls)
allen_10xv3_neu = allen_10xv3_h5ad[allen_10xv3_cl_neu].copy()
allen_10xv3_neu.obs["subclass_id_label"] = allen_10xv3_neu.obs.cl.apply(
    lambda i: cl2subclass[int(i)]
)
allen_10xv3_neu.obs["supertype_id_label"] = allen_10xv3_neu.obs.cl.apply(
    lambda i: cl2supertype[int(i)]
)

# allen_10xv3_neu.write_h5ad(filename = os.path.join(
#     allen_dir, "allen.10xv3.pt.regions.neu.ann.h5ad"
# ))

with open(os.path.join(proj_root, "meta", "AIT21_k8_markers.txt"), "r") as f:
    k8s: List[str] = [g.strip() for g in f.readlines()]

allen_10xv3_neu_k8 = allen_10xv3_neu[:, allen_10xv3_neu.var_names.isin(k8s)].copy()
# Allen use X to store CPM and layer/rawcount as raw count.
allen_10xv3_neu_k8.X = allen_10xv3_neu_k8.layers["rawcount"]
allen_10xv3_neu_k8.layers = None
allen_10xv3_neu_k8.write_h5ad(
    filename=os.path.join(
        allen_dir, "allen.10xv3.pt.regions.neu.imn.ann.k8.rawcount.h5ad"
    )
)

# * perform downsample for neu imn groups
# conda_dir = "/tscc/nfs/home/szu/miniforge3/"
proj_root = gv.pt_projd
# conda_dir = "/home/szu/miniforge3/"
conda_dir = "/tscc/nfs/home/szu/miniforge3/"
conda_bin = os.path.join(conda_dir, "bin/conda")
r_env = "r"
rbin = os.path.join(conda_dir, "envs", r_env, "bin", "Rscript")
sa2_env = "sa2"
tool_dir = os.path.join(proj_root, "extratools", "singlecell")
ann2seurat_script = os.path.join(tool_dir, "ScanpyAnnData2Seurat.R")
allen_dir = os.path.join(proj_root, "data/allen")


def get_ann2seurat_exp(
    annfnm: str, outfnm: str, nds: int = 50, dscol: str = "cl", matGroup: str = "X"
) -> str:
    r = [
        rbin,
        ann2seurat_script,
        f"-f {annfnm}",
        f"-o {outfnm}",
        f"--conda {conda_bin}",
        f"--condaenv {sa2_env}",
        "-d",
        f"--dscol {dscol}",
        f"--nds {nds}",
        f"--matGroup {matGroup}",
    ]
    return " ".join(r)


nds: int = 100
cmd = get_ann2seurat_exp(
    annfnm=os.path.join(
        allen_dir, "allen.10xv3.pt.regions.neu.imn.ann.k8.rawcount.h5ad"
    ),
    outfnm=os.path.join(
        proj_root,
        "data",
        "allen_seurat",
        f"allen.10xv3.pt.regions.neu.imn.k8.cl.ds{nds}.rds",
    ),
    nds=nds,
    dscol="cl",
    matGroup="X",
)
p = Popen(cmd, shell=True)
p.wait()

# * perform downsample for nn imn groups
nn_cls: List[int] = allen_meta[suffix.str.contains("NN|IMN", regex=True)].cl.to_list()
nn_sc: pd.DataFrame = pd.DataFrame(
    [(k, cl2subclass[k]) for k in nn_cls], columns=["cl", "subclass"]
)
nn_sc.insert(column="nt_type_label", value=[cl2ntype[k] for k in nn_cls], loc=2)
allen_10xv3_cl_nn = allen_10xv3_h5ad.obs.cl.apply(lambda x: int(x)).isin(nn_cls)
with open(os.path.join(proj_root, "meta", "AIT21_k8_markers.txt"), "r") as f:
    k8s: List[str] = [g.strip() for g in f.readlines()]
allen_10xv3_nn = allen_10xv3_h5ad[allen_10xv3_cl_nn].copy()
allen_10xv3_nn.obs["subclass_id_label"] = allen_10xv3_nn.obs.cl.apply(
    lambda i: cl2subclass[int(i)]
)
allen_10xv3_nn.obs["supertype_id_label"] = allen_10xv3_nn.obs.cl.apply(
    lambda i: cl2supertype[int(i)]
)

allen_10xv3_nn_k8 = allen_10xv3_nn[:, allen_10xv3_nn.var_names.isin(k8s)].copy()
# Allen use X to store CPM and layer/rawcount as raw count.
allen_10xv3_nn_k8.X = allen_10xv3_nn_k8.layers["rawcount"]
allen_10xv3_nn_k8.layers = None
allen_10xv3_nn_k8.write_h5ad(
    filename=os.path.join(
        allen_dir, "allen.10xv3.pt.regions.nn.imn.ann.k8.rawcount.h5ad"
    )
)

nds = 1000
cmd = get_ann2seurat_exp(
    annfnm=os.path.join(
        allen_dir, "allen.10xv3.pt.regions.nn.imn.ann.k8.rawcount.h5ad"
    ),
    outfnm=os.path.join(
        proj_root,
        "data",
        "allen_seurat",
        f"allen.10xv3.pt.regions.nn.imn.k8.cl.ds{nds}.rds",
    ),
    nds=nds,
    dscol="cl",
    matGroup="X",
)
p = Popen(cmd, shell=True)
p.wait()


# * prepare down-sample direclty on anndata
# 2024-10-04
def downsample_allen(
    ann: ad.AnnData, on: str = "subclass_id_label", nds: int = 200
) -> ad.AnnData:
    use_cells = (
        ann.obs.groupby(on, observed=True)
        .apply(lambda x: x.sample(min(x.shape[0], nds)))
        .index.droplevel(0)
    )
    r = ann[ann.obs_names.isin(use_cells),].copy()
    return r


ann_neu = ad.read_h5ad(
    os.path.join(allend, "allen.10xv3.pt.regions.neu.imn.ann.k8.rawcount.h5ad"),
    backed=None,
)
ann_nn = ad.read_h5ad(
    os.path.join(allend, "allen.10xv3.pt.regions.nn.imn.ann.k8.rawcount.h5ad"),
    backed=None,
)

ann_neuds = downsample_allen(ann_neu, "subclass_id_label", 1000)
ann_neuds.obs.insert(0, "barcode", ann_neuds.obs_names)
ann_neuds.write_h5ad(os.path.join(allend, "allen.10xv3.pt.neu.ds.h5ad"))

ann_nn = ann_nn[~ann_nn.obs.supertype_id_label.str.contains("IMN"),]
ann_nnds = downsample_allen(ann_nn, "subclass_id_label", 10000)
ann_nnds.write_h5ad(os.path.join(allend, "allen.10xv3.pt.nn.ds.h5ad"))

# 2024-11-12
# doubel check the data: no barcode in obs
ann_neu = ad.read_h5ad(
    os.path.join(allend, "allen.10xv3.pt.regions.neu.imn.ann.k8.rawcount.h5ad"),
    backed="r",
)

# 2024-11-18
# 1. put nn and neuron in one place
# 2. downsample
# 3. add class, subclass info if not found
ann_neuds = ad.read_h5ad(
    os.path.join(allend, "allen.10xv3.pt.neu.ds.h5ad"), backed=None)
ann_nnds = ad.read_h5ad(
    os.path.join(allend, "allen.10xv3.pt.nn.ds.h5ad"), backed=None)
annds = ad.concat([ann_nnds, ann_neuds], merge="same")
annds.obs['barcode'] = annds.obs_names
annds.write_h5ad(os.path.join(allend, "allen.10xv3.pt.k8.ds.raw.h5ad"))

# 2025-03-09
# check the number of cells used in co-embedding
# coembed
# 103393
ann_neuds = ad.read_h5ad(
    os.path.join(allend, "allen.10xv3.pt.neu.ds.h5ad"), backed='r')
# 78833
ann_nnds = ad.read_h5ad(
    os.path.join(allend, "allen.10xv3.pt.nn.ds.h5ad"), backed='r')
# 182226
allen_ann = ad.read_h5ad(
    os.path.join(allend, "allen.10xv3.pt.k8.ds.raw.h5ad"), backed='r')

annd = os.path.join(gv.pt_projd, "data")
pt_ann_fnm = os.path.join(annd, "pairedtag_ann", "pt.RNA.ann.k8.ds.raw.h5ad")
snmC_ann_fnm = os.path.join(annd, "snmC_snm3C", "snmC.scaled.gmat.k8.ds.h5ad")
atac_ann_fnm = os.path.join(annd, "snATAC", "snATAC.ptregion.gmat.k8.ds.raw.h5ad")
# 238159
pt_ann = ad.read_h5ad(pt_ann_fnm, backed='r')
# 52094
snmC_ann = ad.read_h5ad(snmC_ann_fnm, backed='r')
# 162569
atac_ann = ad.read_h5ad(atac_ann_fnm, backed='r')
