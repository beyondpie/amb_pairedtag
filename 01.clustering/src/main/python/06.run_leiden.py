import pyprojroot
import os
import sys
from typing import List, Optional
from logging import Logger

import numpy as np
import pandas as pd
from numba.core.errors import NumbaDeprecationWarning
from numba.core.errors import NumbaPendingDeprecationWarning
import warnings
import snapatac2 as sa2

warnings.simplefilter("ignore", category=NumbaDeprecationWarning)
warnings.simplefilter("ignore", category=NumbaPendingDeprecationWarning)

proj_dir = pyprojroot.here()
pydir = os.path.join(proj_dir, "package", "python")
sys.path.insert(0, pydir)
from leiden import run_leiden

# * configs
annfnm: str = snakemake.input["ann"]
out_leiden: str = snakemake.output["leiden"]
out_silht: str = snakemake.output["silht"]
leiden_weight: bool = True if snakemake.params["leiden_weight"] > 0 else False
nleiden: int = snakemake.params["nleiden"]
nsilht: int = snakemake.params["nsilht"]
nsample_silht: int = snakemake.params["nsample_silht"]
threads: int = snakemake.threads
r: float = float(snakemake.wildcards["r"])

# for DEBUG
# work_dir = f"{proj_dir}/01.clustering"
# out_dir = f"{work_dir}/out/test"
# os.makedirs(out_dir, exist_ok=True)
# leiden_weight: bool = True
# threads = 2
# nleiden: int = 3
# nsilht: int = 2
# nsample_silht: int = 50000
# prefix = f"RNA_k8_npc{npc}_k{k_knn}_L1"
# annfnm = os.path.join(out_dir, f"{prefix}_0.h5ad")
# out_leiden = os.path.join(out_dir, f"{prefix}_leiden.csv")
# out_silht = os.path.join(out_dir, f"{prefix}_silht.csv")
# r = 0.1

print(f"Run leiden on reso: {r}.")
print(f"With {nleiden} repeats and {nsilht} silht repeats.")

ann: sa2.AnnData = sa2.read(filename=annfnm, backed="r")
# import pdb
r_leiden, r_silht = run_leiden(
    sa2ann=ann,
    r=r,
    times_r=nleiden,
    times_sil=nsilht,
    logger=None,
    leiden_weight=leiden_weight,
    nsample_sil=nsample_silht,
 )

r_leiden.to_csv(out_leiden, sep=",", header=True, index=False)
r_silht.to_csv(out_silht, sep=",", header=True, index=False)
ann.close()
