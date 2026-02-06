import sys
import os
import pandas as pd
import scanpy as sc
import snapatac2 as sa2

import pyprojroot

proj_root = pyprojroot.here()

i: int = int(sys.argv[1])

dlts: pd.DataFrame = pd.read_csv(
    os.path.join(proj_root, "meta", "est.dlt.rate.20240324.csv"),
    sep=",", header=0
)
dlts.set_index("exp", drop=False, inplace=True)

nfeature = 3000
sa2ann_dir = os.path.join(proj_root, "data", "pairedtag_ann", "snapatac2")
outdlt_dir = os.path.join(proj_root, "00.datapreprocess", "out", "sa2dlt")


def run_sa2dlt(i: int) -> None:
    print(f"{i}")
    r1 = sa2.read(
        filename=os.path.join(sa2ann_dir, f"sa2.AMB_Exp{i}_RNA.bfdlt.h5ad"),
        backed=None
    )
    vg = sc.pp.highly_variable_genes(
        adata=r1, n_top_genes=nfeature, flavor="seurat_v3", inplace=False
    )
    r1.var["selected"] = vg.highly_variable.to_list()
    sa2.pp.scrublet(
        adata=r1,
        features="selected",
        expected_doublet_rate=dlts.loc[f"exp{i}"].dltr,
        inplace=True,
    )
    d = r1.obs[["doublet_probability", "doublet_score"]]
    d.insert(0, "barcode", d.index)
    d.to_csv(
        os.path.join(outdlt_dir, f"sa2.AMB_Exp{i}_RNA.dlt.csv"),
        sep=",",
        index=False,
        header=True,
    )
    return None


if __name__ == "__main__":
    run_sa2dlt(i)
