import os
import sys
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd

import requests
import SimpleITK as sitk
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

import tmpPypkg.globalvar as gv
 
# * meta
abc_atlas = os.path.join(gv.pt_projd, "data", "abc_atlas")
abc_cache = AbcProjectCache.from_cache_dir(abc_atlas)
sects = [f"Zhuang-ABCA-{i}" for i in range(1, 5)]

# * get region information for 3D
# * get MERFISH meta data
# * get MERFISH anndata
# * seperate NN and Neu
# * Region restrictions on Paired-Tag parts
# * perform downsample on subclass level



