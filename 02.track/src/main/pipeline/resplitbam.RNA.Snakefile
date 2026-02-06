envvars:
    "PATH"

import os
import pandas as pd
from typing import List, Dict, Tuple
import pyprojroot

# * configs
nthread = config["nthread"]
debug = config["debug"]

proj_dir = pyprojroot.here()
cell_meta_fnm = config["cellmeta"]
sublib_fnm = config["sublib"]
samtools = config["samtools"]
overwrite = config["overwrite"]
meta_sublib: pd.DataFrame = pd.read_csv(sublib_fnm,
                                        sep = ",", header = 0,
                                        names = None)
sublib_RNAs = meta_sublib['RNA']
sublib_RNA_bams = meta_sublib['RNA_bam']
sublib_ids = meta_sublib['sublibraryid']
d2bam: Dict[str, str] = {
    d: b for d, b in zip(sublib_RNAs, sublib_RNA_bams)}
d2sublib: Dict[str, str] = {
    d: s for d, s in zip(sublib_RNAs, sublib_ids)}

# * set out directories
work_dir = os.path.join(proj_dir, "02.track")
outdir = config['outdir']
flagdir = os.path.join(outdir, "flag")
logdir = os.path.join(outdir, "log")
bamdir = os.path.join(outdir, "cell_bams")
indexdir = os.path.join(outdir, "index_dir")
for d in [outdir, flagdir, logdir, bamdir, indexdir]:
    os.makedirs(d, exist_ok = True)

# sub directories
# column 9
sexs = ['MaleA', 'MaleB', 'FemaleA', 'FemaleB']
# column 4
regions = ['AMY', 'CPU', 'ERC', 'HCa',
           'HCp', 'HYP', 'NAC', 'PFC', 'VTA_SnR']

detailed_bamdirs = [os.path.join(bamdir, f"{r}_{s}_RNA")
                    for r in regions
                    for s in sexs]
for d in detailed_bamdirs:
    os.makedirs(d, exist_ok = True)

# * supplemantary functions
def get_bam_fnm(wildcards) -> str:
    return d2bam[wildcards.d]
def get_sublib_id(wildcards) -> str:
    return d2sublib[wildcards.d]

resplit_sublib_RNAs: List[str] = [
    "ZW253",
    "ZW254",
    "ZW261",
    "ZW262",
    "ZW263",
    "ZW264",
    "ZW265",
    "ZW266",
    "ZW440"
]

rule all:
    input:
        expand("{f}/index_{d}_cell_bams.done", f = flagdir,
               d = resplit_sublib_RNAs),
        expand("{f}/split_{d}_cell_bams.done", f = flagdir,
               d = resplit_sublib_RNAs)

rule index_bam:
    input:
        bam_fnm = get_bam_fnm
    output:
        index_fnm = f"{indexdir}/{{d}}.index",
        flag = touch(f"{flagdir}/index_{{d}}_cell_bams.done")
    threads: 1
    shell:
        """
        {samtools} index {input.bam_fnm} {output.index_fnm}
        """
    

rule split_bam:
    input:
        bam_fnm = get_bam_fnm,
        index_fnm = f"{indexdir}/{{d}}.index",
        cell_meta_fnm = cell_meta_fnm
    params:
        sublib_id = get_sublib_id,
        out_bam_dir = bamdir,
        debug = debug,
        overwrite = overwrite
    threads: nthread
    output:
        flag = touch(f"{flagdir}/split_{{d}}_cell_bams.done")
    log:
        f"{logdir}/split_{{d}}_cell_bams.log"
    threads: nthread
    script:
       f"{work_dir}/src/main/python/02.splitbam_RNA.py" 
