"""
To understand how gimmemotifs works and how can we link our ChromHMM information at this level.
If we cannot directly use this level information, we have to run CellOracle pipeline again for
for the elements we are interested in.
"""

import os
import pandas as pd
import celloracle as co
from celloracle.motif_analysis.tfinfo_core import load_TFinfo

# * meta
projd = "/projects/ps-renlab2/szu/projects/amb_pairedtag/"
cellOracled = os.path.join(projd, "data", "cellOracle")
tfscand = os.path.join(cellOracled, "tfscan")
grnd = os.path.join(cellOracled, "GRN_subclass")

# * functions
def load_tfinfo_from_parquet(f) -> pd.DataFrame:
    import pyarrow.parquet as pq
    r = pq.read_table(f)
    return r.to_pandas()


# * main
testsc = "Astro-TE_NN"
# 1. read tfscan results

tfinfo = load_TFinfo(
    file_path=os.path.join(
        tfscand, f"sa2subclass.{testsc}.pdc.celloracle.tfinfo"))

# 2. check gimmiemotif pipeline
# NOTE: in our snATAC-seq paper, we use 200(default) as background seq size.
# During scan, for each sequence, it will filter out the motifs based on fdr.
# Ideally, each DNA sequence, will have multiple (at least one) motifs after
# fdr.

# our threshold, this should be a further filtering base on
# the motif matching scores.

threshold_score = 10.0
tfinfo.scanned_filtered = tfinfo.scanned_df
tmp = tfinfo.scanned_filtered.groupby(by=["seqname", "motif_id"]).sum()
tmp = tmp[tmp.score >= threshold_score]
tmp = tmp.reset_index()

tfinfo.reset_filtering()
tfinfo.filter_motifs_by_score(threshold=threshold_score)
tfinfo.make_TFinfo_dataframe_and_dictionary(verbose=True)


filter_tfinfo = load_tfinfo_from_parquet(
    os.path.join(tfscand, f"sa2subclass.{testsc}.pdc.baseGRN.df.parquet")
)

# 2. check Cell Oracle results
grn = co.utility.load_hdf5(
    os.path.join(grnd, f"GRN.{testsc}.celloracle.links")
)

# grn.filtered_links
# Out[834]: 
# {'Astro-TE_NN':         source         target  coef_mean  coef_abs             p      -logp
#  79598    Nr2f2        Gm29683   1.172668  1.172668  5.660375e-20  19.247155
#  120623    Mitf         Mgat4c   0.985316  0.985316  1.392738e-16  15.856131
#  1251     Nr2f1  2900052N01Rik   0.807171  0.807171  5.973685e-18  17.223758
#  133166   Tcf12          Nr3c1   0.795614  0.795614  1.192904e-15  14.923394
#  146801   Npas2          Plce1   0.784424  0.784424  1.792522e-15  14.746535
#  ...        ...            ...        ...       ...           ...        ...
#  91421     E2f3          Hmgb2  -0.102844  0.102844  6.330919e-09   8.198533
#  180236  Dmrta2          Smad1   0.102833  0.102833  1.806139e-08   7.743249
#  55771    Npas3           Emp2  -0.102819  0.102819  1.545765e-06   5.810857
#  175312     Id3        Shroom3  -0.102818  0.102818  7.082395e-06   5.149820
#  73293    Stat4           Gfap   0.102816  0.102816  1.393528e-05   4.855884
 
#  [10000 rows x 6 columns]}
