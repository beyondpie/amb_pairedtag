"""
Get single-cell level snmC and snm3C gene activity matrix.

- Here we focus on the major regions in our Paired-Tag data.
- use mCG for NN, IMN and granule cells (DG Glut and CB Granule Glut)
- use mCH for other neurons

- snmC ans snm3C are two different platforms, and the former includes much more
  cells in Paired-Tag region. We firstly focus on snmC.
"""
import os
import re
import xarray as xr
import pandas as pd
import anndata as ad
from ALLCools.mcds import MCDS
import tmpPypkg.globalvar as gv
from tmpPypkg.utils import flatten_list

snmC3Cd = os.path.join(gv.pt_projd, "data/snmC_snm3C", "raw")
pt2CEMBARegion = gv.get_rawCEMBAdissect_in_ambPairedTag()
ptCEMBARegions = set(
    flatten_list([i for i in pt2CEMBARegion.values()]))

# read mcds and m3cds dataset
def load_mcds(fnm: str) -> MCDS:
    r = MCDS.open(
        os.path.join(
            snmC3Cd, fnm),
        obs_dim="cell",
        var_dim="geneslop2k",
        use_obs=None)
    return r
# read metadata for the two modalities
def load_meta(fnm: str) -> pd.DataFrame:
    r = pd.read_csv(
        os.path.join(snmC3Cd, "CEMBA.mC.Metadata.csv"), header=0, sep=",")
    r.set_index('cell', drop=False, inplace=True)
    return r

def select_cell(meta: pd.DataFrame) -> pd.Series:
    return meta[meta.CEMBARegion.isin(ptCEMBARegions)].cell

# from module ALLCools.mcds
def _make_obs_df_var_df(use_data, obs_dim, var_dim):
    obs_df = pd.DataFrame([], index=use_data.get_index(obs_dim).astype(str))
    var_df = pd.DataFrame([], index=use_data.get_index(var_dim).astype(str))
    coord_prefix = re.compile(f"({obs_dim}|{var_dim})_")
    for k, v in use_data.coords.items():
        if k in [obs_dim, var_dim]:
            continue
        try:
            # v.dims should be size 1
            if v.dims[0] == obs_dim:
                series = v.to_pandas()
                # adata.obs_name is str type
                series.index = series.index.astype(str)
                obs_df[coord_prefix.sub("", k)] = series
            elif v.dims[0] == var_dim:
                series = v.to_pandas()
                # adata.var_name is str type
                series.index = series.index.astype(str)
                var_df[coord_prefix.sub("", k)] = series
            else:
                pass
        except IndexError:
            # v.dims is 0, just ignore
            pass
    return obs_df, var_df

def to_ann(da: xr.DataArray, outfnm: str) -> ad.AnnData:
    use_data = da.squeeze()
    obs_df, var_df = _make_obs_df_var_df(use_data, "cell", "geneslop2k")
    adata = ad.AnnData(
        X=use_data.astype('float32').values,
        obs=obs_df,
        var=var_df,
        dtype="float32")
    adata.write_h5ad(outfnm)
    return adata


if __name__ == '__main__':

    mCds = load_mcds("CEMBA.snmC.mcds")
    snmCmeta = load_meta("CEMBA.mC.Metadata.csv")
    # 63792 cells
    snmCells = select_cell(snmCmeta)

    m3Cds = load_mcds("CEMBA.snm3C.mcds")
    snm3Cmeta = load_meta("CEMBA.m3C.Metadata.csv")
    # 11F    1807
    # 9J     1165
    # 9H     1104
    # 11E    1091
    # 8J      883
    # 8E      514
    # So let's focus on snmC firstly
    snm3Cells = select_cell(snm3Cmeta)

    # * select cells and gmat
    raw_cells = mCds.get_index('cell')
    mCH = mCds['geneslop2k_da_frac'].sel(mc_type="CHN")
    mCHpt = mCH.sel(cell=raw_cells.isin(snmCells))

    # load ensmug2genesymbol
    # geneNames = pd.read_csv(
    #     os.path.join(gv.pt_projd, "meta",
    #                 "gencode_GRCm38vM22_id2name.csv"),
    #     sep=",", header=0)
    # geneNames.set_index('gene_id', drop=False, inplace=True)
    # raw_ens = mCds.get_index('geneslop2k')
    # 55487 all in geneNames

    # * to anndata
    adata_mCH = to_ann(
        da=mCHpt,
        outfn=os.path.join(gv.pt_projd, "data/snmC_snm3C",
                        "snmC.mCH.pt.Xonly.h5ad"))

    mCG = mCds['geneslop2k_da_frac'].sel(mc_type="CGN")
    mCGpt = mCG.sel(cell=raw_cells.isin(snmCells))
    adata_mCG = to_ann(
        da=mCGpt,
        outfnm=os.path.join(gv.pt_projd, "data/snmC_snm3C",
                            "snmC.mCG.pt.Xonly.h5ad"))
