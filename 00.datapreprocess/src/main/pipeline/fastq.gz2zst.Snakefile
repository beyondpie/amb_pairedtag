envvars:
    "PATH"

import os

# tscc
#zstd = "/tscc/nfs/home/szu/miniforge3/bin/zstd"

# mediator
zstd = "/home/szu/miniforge3/bin/zstd"

shared_data = "/tscc/projects/ps-renlab2/szu/shared_data"
projd = "/tscc/projects/ps-renlab2/szu/projects/amb_pairedtag"
workd = f"{projd}/00.datapreprocess"
flagd = f"{workd}/out/flag"
logd = f"{workd}/out/log"

with open(f"{projd}/meta/pairedtag.Raw.FASTQgz.R1.csv", 'r') as f:
    lines = [l.strip().split(",") for l in f.readlines()]
    R1FASTQgz_dict = { l[0]: l[1] for l in lines}
    R1FASTQzst_dict = {l[0]: f"{shared_data}/wmb_EP/fastqR/{os.path.basename(l[1])}.zst" for l in lines}

with open(f"{projd}/meta/pairedtag.Raw.FASTQgz.R2.csv", 'r') as f:
    lines = [l.strip().split(",") for l in f.readlines()]
    R2FASTQgz_dict = { l[0]: l[1] for l in lines}
    R2FASTQzst_dict = {l[0]: f"{shared_data}/wmb_EP/fastqR/{os.path.basename(l[1])}.zst" for l in lines}

sublibIds = list(R1FASTQgz_dict.keys())
# test
# sublibIds = list(R1FASTQgz_dict.keys())[0:1]

    
rule all:
    input:
        expand("{f}/{s}_gz2zst_R1.done", f = flagd, s = sublibIds),
        expand("{f}/{s}_gz2zst_R2.done", f = flagd, s = sublibIds)

rule runzstR1:
    input:
        FASTQgz = lambda wildcards: R1FASTQgz_dict[wildcards.s]
    output:
        tag = touch(f"{flagd}/{{s}}_gz2zst_R1.done")
    params:
        FASTQzst = lambda wildcards: R1FASTQzst_dict[wildcards.s]
    resources:
        mem = "20G",
        walltime = "10:00:00",
        queue = "condo",
        mail = "ae",
        email = "debug.pie@gmail.com",
        slurm = "qos=condo account=csd788 nodes=1"
    log:
        f"{logd}/{{s}}_gz2zst_R1.log"
    threads: 20
    shell:
        """
        gzip -d -c {input.FASTQgz} | {zstd} -T{threads} -19 -f -o {params.FASTQzst}
        """
    

rule runzstR2:
    input:
        FASTQgz = lambda wildcards: R2FASTQgz_dict[wildcards.s]
    output:
        tag = touch(f"{flagd}/{{s}}_gz2zst_R2.done")
    params:
        FASTQzst = lambda wildcards: R1FASTQzst_dict[wildcards.s]
    resources:
        mem = "20G",
        walltime = "10:00:00",
        queue = "condo",
        mail = "ae",
        email = "debug.pie@gmail.com",
        slurm = "qos=condo account=csd788 nodes=1"
    log:
        f"{logd}/{{s}}_gz2zst_R2.log"
    threads: 20
    shell:
        """
        gzip -d -c {input.FASTQgz} | {zstd} -T{threads} -19 -f -o {params.FASTQzst}
        """
    

    
    
