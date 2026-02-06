envvars:
    "PATH"
configfile: "config.yaml"
    
import os
from typing import Dict, Any, List
import numpy as np
import pandas as pd
import pyprojroot

proj_dir = pyprojroot.here()
work_dir = f"{proj_dir}/01.clustering"
out_dir = f'{work_dir}/{config["outdir"]}'
fig_dir = f"{out_dir}/figures"
ann_dir = f"{out_dir}/anns"
leiden_dir = f"{out_dir}/leidens"
silht_dir = f"{out_dir}/silhts"
log_dir = f'{out_dir}/log'
flag_dir = f'{out_dir}/flag'

for d in [out_dir, fig_dir,
          ann_dir,
          leiden_dir,
          silht_dir,
          log_dir, flag_dir]:
    os.makedirs(d, exist_ok = True)

# params about leiden
r_start = config['r_start']
r_end = config['r_end']
resos: List[float] = [round(i, 1)
                      for i in np.arange(
                              r_start, r_end + 0.1, 0.1).tolist()]
nleiden: int = config['nleiden']
nsilht: int = config['nsilht']
nsample_silht:int = config['nsample_silht']
nsample_umap: int = config['nsample_umap']
    
cll = config["clustering_level"]
npc = config["npc"]
use_k8: int = config["use_k8"]
nfeature: int = config['nfeature']
k_knn = config["k_knn"]
feat: str = "k8" if use_k8 > 0 else f"vg{nfeature}"
prefix_ann = f"RNA_{feat}_npc{npc}_k{k_knn}_{cll}"

b2g_fnm = config["b2gf"]
ann_h5ad = config["annh5ad"]
allen_k8_fnm = config["allenk8"]

for f in [b2g_fnm, ann_h5ad, allen_k8_fnm]:
    if not os.path.exists(f):
        raise FileNotFoundError(f"No file: {f}.")

b2g:pd.DataFrame = pd.read_csv(b2g_fnm, sep = ',', header = 0)
group2cnt_: Dict[int, int] = b2g[cll].value_counts().to_dict()
# force them to be int
group2cnt = {int(k):int(v) for k, v in group2cnt_.items()}

pre_clusters: List[int] = list(group2cnt.keys())
print(f"{b2g.shape[0]} barcodes with {len(pre_clusters)} cluster(s) in total.")
for k, v in group2cnt.items():
    print(f"group {k}: {v} barcodes.")
nmin_cluster:int = config['nmin_cluster']
print(f"filter pre_cluster: >= {nmin_cluster} cells.")
pre_clusters = [i for i in list(group2cnt.keys()) if group2cnt[i] >= nmin_cluster]
print(f"{len(pre_clusters)} pre_clusters will be considered.")
# add filter clusters
if "filtercl" in config:
    fnm = config["filtercl"]
    filtercl: List[int] = pd.read_csv(fnm, header = None)[0].to_list()
    print(f"{len(filtercl)} found under {fnm}.")
    pre_clusters = [i for i in pre_clusters if i not in filtercl]
    print(f"{len(pre_clusters)} pre_clusters after filtering.")
else:
    print("no filtercl in config.")

rule all:
    input:
        ann = expand("{d}/{f}_c{i}.h5ad",
                     d = ann_dir,
                     f = prefix_ann,
                     i = pre_clusters),
        cellmeta = expand("{d}/{f}_c{i}.cellmeta.csv",
                          d = out_dir,
                          f = prefix_ann,
                          i = pre_clusters),
        leidens = expand("{d}/{f}_c{i}_r{r}.leiden.csv",
                         d = leiden_dir,
                         f = prefix_ann,
                         i = pre_clusters,
                         r = resos),
        silhts = expand("{d}/{f}_c{i}_r{r}.silht.csv",
                        d = silht_dir,
                        f = prefix_ann,
                        i = pre_clusters,
                        r = resos),
        leiden = expand("{d}/{f}_c{i}.leiden.csv",
                        d = out_dir ,
                        f = prefix_ann,
                        i = pre_clusters),
        silht = expand("{d}/{f}_c{i}.silht.csv",
                       d = out_dir,
                       f = prefix_ann,
                       i = pre_clusters),
        umap = expand("{d}/{f}_c{i}.umap.pdf",
               d = fig_dir, f = prefix_ann, i = pre_clusters)

rule sa2_embed:
    input:
        b2g_fnm = b2g_fnm,
        ann_h5ad = ann_h5ad,
        allen_k8 = allen_k8_fnm
    output:
        ann = f"{ann_dir}/{prefix_ann}_c{{i,\d+}}.h5ad",
        cellmeta = f"{out_dir}/{prefix_ann}_c{{i,\d+}}.cellmeta.csv",
        flag = touch(f"{flag_dir}/{prefix_ann}_c{{i,\d+}}.done")
    log:
        f"{log_dir}/{prefix_ann}_c{{i,\d+}}.log"
    conda: "sa2"
    params:
        npc = npc,
        k_knn = k_knn,
        knn_method = "kdtree",
        cll = cll,
        use_k8 = use_k8,
        nfeature = nfeature
    threads: config['ncpu']
    resources:
        time = "15:00:00",
        mem = "200G",
        partition = "platinum",
        slurm="qos=condo account=csd788 nodes=1"
    script:
        f"{work_dir}/src/main/python/05.sa2.embed.py"

rule run_leiden:
    input:
        ann = f"{ann_dir}/{prefix_ann}_c{{i,\d+}}.h5ad",
    output:
        leiden = f"{leiden_dir}/{prefix_ann}_c{{i}}_r{{r}}.leiden.csv",
        silht = f"{silht_dir}/{prefix_ann}_c{{i}}_r{{r}}.silht.csv"
    conda: "sa2"
    threads: 1
    resources:
        time = "07:00:00",
        mem = "15G",
        partition = "platinum",
        slurm="qos=condo account=csd788 nodes=1"
    params:
        leiden_weight = 1,
        nleiden = nleiden,
        nsilht = nsilht,
        nsample_silht = nsample_silht
    script:
        f"{work_dir}/src/main/python/06.run_leiden.py"

        
rule sum_leiden:
    input:
        leidens = expand("{o}/{p}_c{{i}}_r{r}.leiden.csv",
                        o = leiden_dir, p = prefix_ann,
                        r = resos),
        silhts = expand("{o}/{p}_c{{i}}_r{r}.silht.csv",
                        o = silht_dir, p = prefix_ann,
                        r = resos),
        cellmeta = f"{out_dir}/{prefix_ann}_c{{i}}.cellmeta.csv"
    output:
        leiden = f"{out_dir}/{prefix_ann}_c{{i}}.leiden.csv",
        silht = f"{out_dir}/{prefix_ann}_c{{i}}.silht.csv",
        cellmeta = f"{out_dir}/{prefix_ann}_c{{i}}.cellmeta.leiden.csv",
        umap = f"{fig_dir}/{prefix_ann}_c{{i}}.umap.pdf"
    conda: "r"
    params:
        nsample_umap = nsample_umap
    threads: 1
    resources:
        time = "02:00:00",
        mem = "50G",
        partition = "platinum",
        slurm="qos=condo account=csd788 nodes=1"
    script:
        f"{work_dir}/src/main/R/07.sum_leiden.R"
