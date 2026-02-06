import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from pybedtools import BedTool
from sklearn.preprocessing import RobustScaler, StandardScaler
import warnings


def scanpy_PCA_plus(
        a: ad.AnnData, n_comps: int,
        weight_by_var: bool = True, **kwargs) -> None:
    """
    Ref:
    https://github.com/satijalab/seurat/blob/1549dcb3075eaeac01c925c4b4bb73c73450fc50/R/dimensional_reduction.R#L897
    """
    sc.pp.pca(a, n_comps=n_comps, **kwargs)
    key_obsm = "X_pca"
    key_uns = "pca"
    if not weight_by_var:
        print("Normalize PCA by diving the singluar values.")
        sdev = np.sqrt(a.uns[key_uns]['variance'])
        singular_values = sdev * np.sqrt(a.shape[0] - 1)
        old_X_pca = a.obsm[key_obsm]
        new_X_pca = old_X_pca / singular_values
        a.obsm[key_obsm] = new_X_pca
        a.obsm[f"scanpy_{key_obsm}"] = old_X_pca
        a.uns[key_uns]['singular_values'] = singular_values
    else:
        print("Keep the scanpy default PCA, i.e., weighted by variance of PCs.")
    return None

def normalize_data(adata: ad.AnnData,
                           target_sum: float = 1e6) -> None:
    print("Save raw counts into layer.")
    adata.layers['raw_counts'] = adata.X
    print("Perform logCPM normalization for raw_counts.")
    sc.pp.normalize_total(adata,
                          target_sum=target_sum,
                          exclude_highly_expressed=False,
                          max_fraction=0.05,
                          key_added=None,
                          inplace=True,
                          layer="raw_counts",
                          copy=False)
    sc.pp.log1p(adata, base=np.e)
    print("Scaling the logCPM.")
    sc.pp.scale(adata, zero_center=True)
    return None

def downsample_ann(ann: ad.AnnData,
                   on: str = "subclass_id_label",
                   nds: int = 200) -> ad.AnnData:
    use_cells = ann.obs.groupby(on, observed=True).apply(
        lambda x: x.sample(min(x.shape[0], nds))).index.droplevel(0)
    r = ann[ann.obs_names.isin(use_cells), ].copy()
    return r


def remove_black_list_region(adata, black_list_path, region_axis=1, f=0.2):
    """
    Remove regions overlap (bedtools intersect -f {f}) with regions in the black_list_path.

    Parameters
    ----------
    adata
        AnnData object
    black_list_path
        Path to the black list bed file
    region_axis
        Axis of regions. 0 for adata.obs, 1 for adata.var
    f
        Fraction of overlap when calling bedtools intersect
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if region_axis == 1:
            feature_bed_df = adata.var[["chrom", "start", "end"]]
        elif region_axis == 0:
            feature_bed_df = adata.obs[["chrom", "start", "end"]]
        else:
            raise ValueError("region_axis should be 0 or 1.")
        feature_bed = BedTool.from_dataframe(feature_bed_df)
        black_list_bed = BedTool(black_list_path)
        black_feature = feature_bed.intersect(black_list_bed, f=f, wa=True)
        try:
            black_feature_index = (
                black_feature.to_dataframe().set_index(["chrom", "start", "end"]).index
            )
            black_feature_id = pd.Index(
                feature_bed_df.reset_index()
                .set_index(["chrom", "start", "end"])
                .loc[black_feature_index][feature_bed_df.index.name]
            )
            print(
                f"{black_feature_id.size} regions removed due to overlapping"
                f" (bedtools intersect -f {f}) with black list regions."
            )
            if region_axis == 1:
                adata._inplace_subset_var(~adata.var_names.isin(black_feature_id))
            else:
                adata._inplace_subset_obs(~adata.obs_names.isin(black_feature_id))
        except pd.errors.EmptyDataError:
            print("no overlap with black list")
            pass
    return


def remove_chromosomes(
    adata, exclude_chromosomes=None,
    include_chromosomes=None, chrom_col="chrom"):
    """Remove chromosomes from adata.var."""
    judge = None
    if exclude_chromosomes is not None:
        not_to_exclude = ~adata.var[chrom_col].isin(exclude_chromosomes)
        judge = not_to_exclude
    if include_chromosomes is not None:
        include = adata.var[chrom_col].isin(include_chromosomes)
        if judge is None:
            judge = include
        else:
            judge &= include

    if judge is not None:
        adata._inplace_subset_var(judge)
        print(f"{adata.shape[1]} regions remained.")
    return

def log_scale(
    adata, method="standard", with_mean=False, with_std=True, max_value=10, scaler=None
):
    """
    Perform log transform and then scale the cell-by-feature matrix.

    Parameters
    ----------
    adata
        adata with normalized, unscaled cell-by-feature matrix
    method
        the type of scaler to use:
        'standard' for :class:`sklearn.preprocessing.StandardScaler`;
        'robust' for :class:`sklearn.preprocessing.RobustScaler`.
    with_mean
        Whether scale with mean center
    with_std
        Whether scale the std
    max_value
        Whether clip large values after scale
    scaler
        A fitted sklearn scaler, if provided, will only use it to transform the adata.

    Returns
    -------
    adata.X is scaled in place, the fitted scaler object will be return if the `scaler` parameter is None.
    """
    # log transform
    if "log" in adata.uns and adata.uns["log"]:
        # already log transformed
        print("adata.X is already log transformed, skip log step.")
    else:
        adata.X = np.log(adata.X)
        adata.uns["log"] = True

    if scaler is not None:
        print("Using user provided scaler.")
        if isinstance(scaler, str):
            import joblib

            scaler = joblib.load(scaler)
        user_provide = True
    else:
        if method == "robust":
            scaler = RobustScaler(with_centering=with_mean, with_scaling=with_std)
        else:
            scaler = StandardScaler(with_mean=with_mean, with_std=with_std)
        user_provide = False

    # transform data
    if user_provide:
        adata.X = scaler.transform(adata.X)
    else:
        adata.X = scaler.fit_transform(adata.X)

    # clip large values
    if max_value is not None:
        adata.X[adata.X > max_value] = max_value
        adata.X[adata.X < -max_value] = -max_value

    # return scaler or not
    if user_provide:
        return
    else:
        return scaler

def neg_snmC(adata: ad.AnnData) -> None:
    if "neg_snmC" in adata.uns and adata.uns["neg_snmC"]:
        print("adata.X has taken negative mC level.")
        return None
    adata.X *= -1
    adata.uns["neg_snmC"] = True
    return None
