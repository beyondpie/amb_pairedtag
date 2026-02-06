import argparse
import datetime
import warnings
import pandas as pd
import polars as pl
import snapatac2 as sa2

warnings.filterwarnings("ignore")
def getDE(ann: sa2.AnnData,
          g1: str,
          outfnm: str,
          groupby: str = "annot.sc",
          min_pct: float = 0.005,
          min_log2fc: float = 0.1,
          ds_fold: float = 2.0) -> pd.DataFrame:
    """
    NOTE: ann X is raw count. 
    """
    print(f"start at {datetime.datetime.now()}.")
    n_g1 = sum(ann.obs[groupby] == g1)
    n_g2 = sum(ann.obs[groupby] != g1)
    g1_barcodes = ann.obs_names[ann.obs[groupby] == g1].to_series()
    g2_barcodes = (
        ann.obs_names[ann.obs[groupby] != g1]
        .to_series()
        .sample(n=min(n_g2, int(ds_fold * n_g1)), replace=False)
    )
    diff_peaks: pl.DataFrame = sa2.tl.diff_test(
        data=ann,
        cell_group1=g1_barcodes,
        cell_group2=g2_barcodes,
        direction="positive",
        min_pct=min_pct,
        min_log_fc=min_log2fc)
    print(f"out result to {outfnm}.")
    diff_peaks.write_csv(outfnm, include_header=True, separator="\t")
    return diff_peaks

# * debug
# import os
# ann = sa2.read(
#     os.path.join("/projects/ps-renlab2/szu/projects/amb_pairedtag",
#                  "12.DE", "out", "ann.H3K9me3.pmat.h5ad"),
#     backed=None
# )

# g1 = "001 CLA-EPd-CTX Car3 Glut"
# groupby = "annot.sc"
# min_pct = 0.005
# min_log2fc = 0.5
# ds_fold = 2.0
# outfnm = "haha.tsv"


# * main
parser = argparse.ArgumentParser(
    prog="DE",
    description="Differential peak analysis by SnapATAC2 logistic model",
)
parser.add_argument("--annfnm", type=str, help="scanpy AnnData with raw count.")
parser.add_argument("--g1", type=str, help="subclass for differential test")
parser.add_argument(
    "--groupby", type=str, default="annot.sc", help="column in ann for group"
)
parser.add_argument(
    "--minpct", type=float, default=0.005, help="min percentage for filtering peaks"
)
parser.add_argument(
    "--minlog2fc",
    type=float,
    default=0.1,
    help="min log2 fold change for filtering peaks",
)
parser.add_argument(
    "--dsfold",
    type=float,
    default=5.0,
    help="n fold for downsampling to prepare background",
)
parser.add_argument("--outfnm", type=str, help="output of results")
args = parser.parse_args()
ann = sa2.read(args.annfnm, backed=None)
r = getDE(
    ann,
    g1=args.g1,
    outfnm=args.outfnm,
    groupby=args.groupby,
    min_pct=args.minpct,
    min_log2fc=args.minlog2fc,
    ds_fold=args.dsfold
)
