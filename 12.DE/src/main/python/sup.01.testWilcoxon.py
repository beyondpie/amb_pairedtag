import os
import datetime

import numpy as np
import pandas as pd

import scanpy as sc
import anndata as ad

# * meta
projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag"
sa2h5adir = os.path.join("/projects/ps-renlab2/zhw063",
                         "99.MouseBrainPairedTag",
                         "20241202_for_Songpeng")
outd = os.path.join(projd, "12.DE", "out")

# * perform normalization
ann_h27ac = ad.read_h5ad(
    os.path.join(outd, "ann.H3K27ac.pmat.bedtool.h5ad"), backed=None)

# * DE analysis
def getPseudoBulkCPM(ann: ad.AnnData, target_sum=1e6) -> np.ndarray:
    """Get pseudobulk CPM for one group
    """
    r = np.sum(ann.X, axis=0)
    n_total = np.sum(r)
    return np.asarray((r / n_total) * target_sum).squeeze()


def getDEbyWilcoxon(ann: ad.AnnData,
                    g1: str,
                    groupby: str = "annot.sc",
                    min_pct: float = 0.005,
                    min_log2fc: float = 0.1,
                    ds_fold: float = 5.0,
                    outfnm: str = "") -> pd.DataFrame:
    """
    NOTE: ann'x is treated as raw count.
    """
    print(f"start at {datetime.datetime.now()}.")
    # perform downsample
    g1_barcodes = ann.obs_names[ann.obs[groupby] == g1].to_series()
    n_g1 = sum(ann.obs[groupby] == g1)
    n_g2 = sum(ann.obs[groupby] != g1)
    # downsample g2 cells
    g2_barcodes = ann.obs_names[
        ann.obs[groupby] != g1].to_series().sample(
            n=min(n_g2, int(ds_fold * n_g1)), replace=False)
    # filter features by min_pct and log2fc
    g1_pct = pd.Series(
        data=np.asarray(
            np.sum(ann[g1_barcodes].X > 0, axis=0)).squeeze() / n_g1,
        index=ann.var_names)
    fea_by_pct = ann.var_names[g1_pct >= min_pct]
    print(
        f"{(fea_by_pct).size} features left after min_pct as {min_pct}.")
    g2_pct = pd.Series(
        data=np.asarray(
            np.sum(ann[g2_barcodes].X > 0, axis=0)).squeeze() / n_g2,
        index=ann.var_names)
    # filter features by min_log2fc
    g1_cpm = getPseudoBulkCPM(ann[g1_barcodes])
    g2_cpm = getPseudoBulkCPM(ann[g2_barcodes])
    log2_fc = pd.Series(np.log2((g1_cpm + 1e-9) / (g2_cpm + 1e-9)),
                        index=ann.var_names)
    fea_by_fc = ann.var_names[log2_fc >= min_log2fc]
    fea = fea_by_pct.intersection(fea_by_fc)
    print(
      f"{fea.size} features left after min_log2fc as {min_log2fc}."
    )
    ann1 = ann[pd.concat([g1_barcodes, g2_barcodes]), ].copy()
    print("Normalize ann: CP10K per cell.")
    sc.pp.normalize_total(ann1, target_sum=1e4, inplace=True)
    # run wilcoxon
    print(f"Run Wilcoxon test on group: {g1}.")
    r = sc.tl.rank_genes_groups(
        ann1[:, fea],
        groupby=groupby,
        groups=[g1],
        reference="rest",
        pts=True,
        copy=True,
        method="wilcoxon"
    )
    print(f"end at {datetime.datetime.now()}.")
    # transform result
    rde = r.uns['rank_genes_groups']
    peaks = np.array([i[0] for i in rde['names']])
    pval = np.array([i[0] for i in rde['pvals']])
    padj = np.array([i[0] for i in rde['pvals_adj']])
    padj = np.array([i[0] for i in rde['pvals_adj']])
    pct = g1_pct[peaks]
    rest_pct = g2_pct[peaks]
    # use my foldchange
    logfd = log2_fc[pd.Series(peaks)]
    rr = pd.DataFrame(
        data={
            'peak': peaks,
            'log2fd': logfd,
            'pct': pct,
            'restpct': rest_pct,
            'pval': pval,
            'padj': padj
        },
        index=peaks
    )
    rr = rr.sort_values(by='padj', ascending=True)
    if len(outfnm) > 1:
        print(f"out result to {outfnm}.")
        rr.to_csv(
            outfnm, header=True, index=False, sep="\t")
    return rr


rr = getDEbyWilcoxon(
    ann=ann_h27ac,
    g1="319 Astro-TE NN",
    groupby="annot.sc",
    min_pct=0.005,
    min_log2fc=0.1,
    ds_fold=5.0,
    outfnm=os.path.join(outd, "DE.319AscTENN.pct-0.005_logfd-0.1.tsv")
)
