import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import snapatac2 as snap

brain_h3k27ac = sc.read_h5ad("/tscc/projects/ps-renlab2/zhw063/99.MouseBrainPairedTag/2024Mar_new_mtxs/20240405.brain.all.H3K27ac.h5ad")
snap.pp.select_features(brain_h3k27ac, n_features=250000, 
                             filter_lower_quantile=0.05, 
                             filter_upper_quantile=0.05, 
                             whitelist=None, 
                             blacklist='/tscc/projects/ps-renlab2/zhw063/99.MouseBrainPairedTag/seurat_v5_objects/mm10.blacklist.bed.gz', 
                             max_iter=1, 
                             inplace=True, 
                             n_jobs=16, 
                             verbose=True)

criteria = brain_h3k27ac.var.index.str.startswith('chrX') | brain_h3k27ac.var.index.str.startswith('chrY') | brain_h3k27ac.var.index.str.startswith('chrM')
brain_h3k27ac.var.loc[criteria, 'selected'] = False

## Everything till here is fine

a = snap.tl.spectral(brain_h3k27ac, features = 'selected', inplace = True)

