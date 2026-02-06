import os
import sys
import numpy as np
import pandas as pd
import polars as pl
import snapatac2 as sa2
import anndata as ad

import tmpPypkg.globalvar as gv
from tmpPypkg.utils import reimport
reimport("tmpPypkg")
reimport("tmpPypkg.globalvar")

# * load snapatac2 object with fragment (no need to load into memory)
atac_sa2annset = sa2.read_dataset(
  filename=gv.atacsa2_anndataset_fnm, mode='r+')

atac_sa2annset.obs['barcode'] = atac_sa2annset.obs_names
atac_sa2annset.close()

# * output the bed.zstd for each singlecell
# the codes below works, but we need to focus on our paired-tag regions
# def chunk_list(lst, size: int = 1000):
#   return [lst[i:i + size] for i in range(0, len(lst), size)]

# atac_sa2annset = sa2.read_dataset(
#   filename=gv.atacsa2_anndataset_fnm, mode='r')

# barcode_chunk = chunk_list(atac_sa2annset.obs_names, size=500)

# for i, bc in enumerate(barcode_chunk):
#   print(f"Running barcode group: {i}")
#   sa2.ex.export_fragments(
#     adata=atac_sa2annset,
#     groupby='barcode',
#     selections=bc,
#     out_dir=os.path.join(gv.pt_projd, "data", "snATAC", "scbedzst")
#   )
#   print(f"Finish running barcode group: {i}")

# * subset to dissections in Paired-Tag
region2disect = gv.get_snATACdissect_in_ambPairedTag()
snATAC_dissect = list(
    set([i for _, v in region2disect.items() for i in v]))

atac_sa2annset_all = sa2.read_dataset(
    filename=gv.atacsa2_anndataset_fnm, mode='r+')

# this step takes about 2 hour and >300G (or 600?) memory
atac_sa2annset_all.obsm['fragment_paired'] = atac_sa2annset_all.adatas.obsm[
    'fragment_paired']

atac_sa2annset_all.to_adata(
    file=os.path.join(gv.pt_projd, "data/snATAC", "all.sa2v26.ann.h5ad"))
atac_sa2annset_all.close()

# focus on ptregion
atac_sa2ann_all = sa2.read(
    os.path.join(gv.pt_projd, "data/snATAC", "all.sa2v26.ann.h5ad"), backed='r')

dissect = pd.Series([
    s.split("_")[1] for s in atac_sa2ann_all.obs["sample"]])
obs_index = dissect.isin(snATAC_dissect)

# takes about 0.5 hours and < 100G memory
# scanpy AnnData
atac_sa2ann_pt = atac_sa2ann_all.subset(
    obs_indices=obs_index, inplace=False)

atac_sa2ann_pt.write(
    os.path.join(gv.pt_projd, "data/snATAC", "atac.ann.ptregion.h5ad"))

# now load as snapatac2 AnnData and add obs information
atac_sa2ann_pt = sa2.read(
    os.path.join(gv.pt_projd, "data/snATAC", "atac.ann.ptregion.h5ad"),
    backed='r+')

annpmat = sa2.read(
    os.path.join(gv.pt_projd,
                 "data/snATAC", "all.sa2v26.ann.pmat.nofrag.h5ad"),
    backed="r")

obs = pd.DataFrame.from_dict(
    {
        'subclass': annpmat.obs['subclass'],
        'cluster': annpmat.obs['cluster'],
        'ptgroup': annpmat.obs['ptgroup']
     }
)
obs.insert(0, 'barcode', annpmat.obs_names)
obs.set_index('barcode', drop=False, inplace=True)

atac_meta_pt = obs.loc[atac_sa2ann_pt.obs_names]
atac_sa2ann_pt.obs['subclass'] = atac_meta_pt.subclass
atac_sa2ann_pt.obs['cluster'] = atac_meta_pt.cluster
atac_sa2ann_pt.obs['ptgroup'] = atac_meta_pt.ptgroup
atac_sa2ann_pt.close()

# * perform subclass bed gz
ann = sa2.read(
    os.path.join(gv.pt_projd, "data/snATAC", "atac.ann.ptregion.h5ad"),
    backed='r')

# about half an hour
sa2.export.export_fragments(
    ann, groupby='subclass',
    out_dir=os.path.join(
        gv.pt_projd, "data", "snATAC", "subclass_bed_zst"))
ann.close()
