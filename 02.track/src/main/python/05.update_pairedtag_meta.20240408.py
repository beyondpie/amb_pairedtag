import os
import pandas as pd
import pyprojroot

proj_dir = pyprojroot.here()
out_cellmeta_fnm = os.path.join(
    proj_dir, "meta", "pairedtag.meta.csv",
)

from_cellmeta_fnm = os.path.join(
    proj_dir, "meta", "pt.barcode.meta.with.dlt.20240330.csv"
)

cell_meta = pd.read_csv(
    from_cellmeta_fnm, sep = ",", header = 0,
    index_col = False
)

r1 = cell_meta.drop(columns = ['exp_x'])
r2 = r1.rename(columns = {'exp_y': 'exp'})
r2.to_csv(out_cellmeta_fnm, sep = ",", header = True, index = False)
